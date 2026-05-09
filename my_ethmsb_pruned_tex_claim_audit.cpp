#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string tex_path;
    std::string out_path;
};

struct Issue {
    std::string id;
    int line = 0;
    std::string text;
    std::string problem;
    std::string rewrite;
};

Options parse(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--tex")
            opt.tex_path = need("--tex");
        else if (a == "--out")
            opt.out_path = need("--out");
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    return opt;
}

std::string default_tex_path()
{
    for (const char* p :
         {"Chapter_3_revised.tex", "Chapter_3_revised (1).tex",
          "/home/wrn/paper tex/Chapter_3_revised.tex"}) {
        std::ifstream in(p);
        if (in) return p;
    }
    throw std::runtime_error("Chapter_3_revised.tex not found");
}

bool has(const std::string& s, const std::string& needle)
{
    return s.find(needle) != std::string::npos;
}

std::string csv_quote(std::string s)
{
    std::string out = "\"";
    for (const char c : s) {
        if (c == '"') out += "\"\"";
        else if (c == '\t')
            out += ' ';
        else
            out += c;
    }
    out += '"';
    return out;
}

void add(std::vector<Issue>& issues, const std::string& id, const int line,
         const std::string& text, const std::string& problem,
         const std::string& rewrite)
{
    issues.push_back({id, line, text, problem, rewrite});
}

std::vector<Issue> audit(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open tex: " + path);
    std::vector<Issue> issues;
    std::string line;
    int line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (has(line, "正确性") && !has(line, "失败率") &&
            !has(line, "概率") && !has(line, "经验"))
            add(issues, "correctness_without_failure_probability", line_no,
                line,
                "single-round variant mentions correctness without empirical "
                "failure probability",
                "将正确性改写为经验失败率约束，并给出随机与边界测试的 Wilson 置信区间。");
        if ((has(line, "稳定") || has(line, "可靠")) && has(line, "提取"))
            add(issues, "stable_reliable_single_round_bitextract", line_no,
                line,
                "single-round BitExtract is boundary-sensitive",
                "单轮短周期 LUT 在边界附近可能失败，本文以失败率实测和参数选择控制该风险。");
        if ((has(line, "p-5") || has(line, "p - 5") || has(line, "p−5")) &&
            !has(line, "p=33") && !has(line, "不支持") &&
            !has(line, "失败"))
            add(issues, "unsafe_pminus5_default", line_no, line,
                "default p-5 needs p=33 unsupported/unsafe caveat",
                "默认候选由 sweep 选择；p=33 下 k=p-5 在当前参数中不可作为 correctness-supported default。");
        if (has(line, "100%") && has(line, "正确"))
            add(issues, "zero_failure_claim", line_no, line,
                "zero-failure correctness claim is not allowed for this branch",
                "改为报告 empirical failure rate，并明确非 worst-case correctness。");
        if (has(line, "可证明安全") &&
            (has(line, "正确") || has(line, "比较")))
            add(issues, "security_vs_correctness_conflation", line_no, line,
                "semantic preservation of CMUX pruning must not be conflated "
                "with comparator correctness",
                "CMUX skip 是语义保持的：若 Rot_a(TV)=TV，则跳过该 CMUX 不改变该 PBS 的明文函数。");
        if ((has(line, "4-32") || has(line, "4 到 32") ||
             has(line, "4至32")) &&
            has(line, "正确"))
            add(issues, "missing_boundary_ci_for_4_32", line_no, line,
                "4-32 bit correctness statement needs boundary tests and CI",
                "周期剪枝比较在各位宽下报告边界失败率、随机失败率与 Wilson 95% 上界。");
        if (has(line, "99.89%"))
            add(issues, "empirical_accuracy_not_worst_case", line_no, line,
                "99.89% must be labeled empirical random accuracy",
                "表中 99.89% 是 empirical random accuracy，不是 worst-case correctness。");
        if (has(line, "TruncRepeat") || has(line, "MetaPBS") ||
            has(line, "迭代式纠偏") || has(line, "逐步吸收残余偏移"))
            add(issues, "metapbs_truncrepeat_out_of_scope", line_no, line,
                "current implementation scope excludes MetaPBS/TruncRepeat",
                "本文实现限定为单轮周期剪枝 BitExtract；多轮纠偏仅可作为 future work，不能作为当前结果依据。");
        if (has(line, "高精度") && (has(line, "MSB") || has(line, "比较")))
            add(issues, "high_precision_needs_probability_context", line_no,
                line,
                "high-precision wording needs probabilistic context for "
                "single-round variant",
                "在 BitExtract 与 mask conversion 正确的条件下，gap offset 将 final PBS 阈值容忍度增加到 (w_k+1)Δ/2。");
    }
    return issues;
}

void emit(std::ostream& out, const std::vector<Issue>& issues)
{
    out << "claim_issue_id,line_number,current_text,problem,"
           "suggested_rewrite\n";
    for (const Issue& i : issues) {
        out << i.id << ',' << i.line << ',' << csv_quote(i.text) << ','
            << csv_quote(i.problem) << ',' << csv_quote(i.rewrite) << "\n";
    }
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        Options opt = parse(argc, argv);
        if (opt.tex_path.empty()) opt.tex_path = default_tex_path();
        const std::vector<Issue> issues = audit(opt.tex_path);
        std::ofstream file;
        std::ostream* out = &std::cout;
        if (!opt.out_path.empty()) {
            file.open(opt.out_path);
            if (!file) throw std::runtime_error("cannot open audit output");
            out = &file;
        }
        emit(*out, issues);
        std::cerr << "PRUNED_TEX_CLAIM_AUDIT path=" << opt.tex_path
                  << " issues=" << issues.size() << "\n";
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: tex_claim_audit\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr
            << "  ./build-ethmsb-release/my_ethmsb_pruned_tex_claim_audit\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
