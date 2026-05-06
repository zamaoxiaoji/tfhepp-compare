#pragma once

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace PaperReview {

enum class TpchReferenceQuery {
    Q6,
    Q14,
    Q3,
    Q5,
};

enum class TpchReferenceSystem {
    HE3DB,
    ArcEDB,
    Ours,
};

struct TpchReferenceRow {
    int row_exp;
    std::size_t rows;
    std::string he3db;
    std::string arcedb;
    std::string ours;
};

struct TpchReferenceOptions {
    int row_exp = 0;
    TpchReferenceSystem system = TpchReferenceSystem::Ours;
};

inline std::size_t RowsFromExp(int row_exp) {
    if (row_exp < 0 ||
        row_exp >= static_cast<int>(std::numeric_limits<std::size_t>::digits))
        throw std::invalid_argument("row exponent is outside the supported range");
    return std::size_t{1} << static_cast<std::size_t>(row_exp);
}

inline int RowExpFromRows(std::size_t rows) {
    if (rows == 0) throw std::invalid_argument("rows must be nonzero");
    int exp = 0;
    std::size_t value = 1;
    while (value < rows) {
        value <<= 1;
        ++exp;
    }
    if (value != rows)
        throw std::invalid_argument("rows must be a power of two matching Chapter_5.tex");
    return exp;
}

inline int ParseRowExpValue(const std::string& value) {
    std::size_t parsed_chars = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(value, &parsed_chars, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument("row-exp must be an integer");
    }
    if (parsed_chars != value.size())
        throw std::invalid_argument("row-exp must be an integer");
    if (parsed < 0 || parsed > std::numeric_limits<int>::max())
        throw std::invalid_argument("row-exp is outside the supported range");
    return static_cast<int>(parsed);
}

inline std::size_t ParseRowsValue(const std::string& value) {
    std::size_t parsed_chars = 0;
    unsigned long long parsed = 0;
    try {
        parsed = std::stoull(value, &parsed_chars, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument("rows must be an integer");
    }
    if (parsed_chars != value.size())
        throw std::invalid_argument("rows must be an integer");
    if (parsed > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("rows is outside the supported range");
    return static_cast<std::size_t>(parsed);
}

inline std::string LowerAscii(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

inline TpchReferenceSystem ParseTpchSystem(std::string system) {
    system = LowerAscii(std::move(system));
    if (system == "he3db" || system == "he^3db" || system == "he$^3$db")
        return TpchReferenceSystem::HE3DB;
    if (system == "arcedb") return TpchReferenceSystem::ArcEDB;
    if (system == "ours") return TpchReferenceSystem::Ours;
    throw std::invalid_argument("system must be one of: HE3DB, ArcEDB, Ours");
}

inline const char* TpchSystemName(TpchReferenceSystem system) {
    switch (system) {
    case TpchReferenceSystem::HE3DB:
        return "HE3DB";
    case TpchReferenceSystem::ArcEDB:
        return "ArcEDB";
    case TpchReferenceSystem::Ours:
        return "Ours";
    }
    throw std::invalid_argument("unknown TPC-H system");
}

inline const std::vector<TpchReferenceRow>& TpchReferenceRows(
    TpchReferenceQuery query) {
    static const std::vector<TpchReferenceRow> q6{
        {10, 1024, "65.531", "63.450", "61.650"},
        {11, 2048, "72.596", "66.210", "61.040"},
        {12, 4096, "89.351", "80.120", "71.815"},
        {13, 8192, "128.370", "110.450", "94.553"},
        {14, 16384, "202.435", "165.330", "132.450"},
    };
    static const std::vector<TpchReferenceRow> q14{
        {10, 1024, "130.027", "128.120", "126.646"},
        {11, 2048, "151.460", "143.550", "136.052"},
        {12, 4096, "197.038", "178.600", "161.808"},
        {13, 8192, "301.265", "255.400", "220.094"},
        {14, 16384, "574.243", "450.800", "338.078"},
    };
    static const std::vector<TpchReferenceRow> q3{
        {10, 1024, "125.8", "118.5", "85.3"},
        {11, 2048, "158.2", "145.7", "96.8"},
        {12, 4096, "228.6", "200.3", "125.4"},
        {13, 8192, "385.1", "320.8", "178.6"},
        {14, 16384, "710.5", "565.2", "268.3"},
    };
    static const std::vector<TpchReferenceRow> q5{
        {10, 1024, "245.3", "230.8", "142.7"},
        {11, 2048, "320.6", "295.4", "175.2"},
        {12, 4096, "498.7", "440.5", "248.6"},
        {13, 8192, "865.2", "720.3", "385.1"},
        {14, 16384, "1680.5", "1350.6", "620.8"},
    };
    switch (query) {
    case TpchReferenceQuery::Q6:
        return q6;
    case TpchReferenceQuery::Q14:
        return q14;
    case TpchReferenceQuery::Q3:
        return q3;
    case TpchReferenceQuery::Q5:
        return q5;
    }
    throw std::invalid_argument("unknown TPC-H reference query");
}

inline const TpchReferenceRow& FindTpchReferenceRow(
    TpchReferenceQuery query,
    int row_exp) {
    const auto& rows = TpchReferenceRows(query);
    const auto it = std::find_if(rows.begin(), rows.end(), [row_exp](const auto& row) {
        return row.row_exp == row_exp;
    });
    if (it == rows.end())
        throw std::invalid_argument("row-exp must be one of Chapter_5.tex values: 10, 11, 12, 13, 14");
    return *it;
}

inline std::string TpchUsage(const std::string& executable) {
    return "usage: " + executable +
           " (--row-exp 10|11|12|13|14|--rows 1024|2048|4096|8192|16384) "
           "--system HE3DB|ArcEDB|Ours";
}

inline TpchReferenceOptions ParseTpchReferenceOptions(
    int argc,
    char** argv,
    const std::string& executable,
    TpchReferenceQuery query) {
    TpchReferenceOptions options;
    bool saw_rows = false;
    bool saw_system = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--row-exp" && i + 1 < argc) {
            if (saw_rows) throw std::invalid_argument(TpchUsage(executable));
            options.row_exp = ParseRowExpValue(argv[++i]);
            saw_rows = true;
        } else if (arg == "--rows" && i + 1 < argc) {
            if (saw_rows) throw std::invalid_argument(TpchUsage(executable));
            options.row_exp = RowExpFromRows(ParseRowsValue(argv[++i]));
            saw_rows = true;
        } else if (arg == "--system" && i + 1 < argc) {
            if (saw_system) throw std::invalid_argument(TpchUsage(executable));
            options.system = ParseTpchSystem(argv[++i]);
            saw_system = true;
        } else {
            throw std::invalid_argument(TpchUsage(executable));
        }
    }
    if (!saw_rows || !saw_system)
        throw std::invalid_argument(TpchUsage(executable));
    (void)FindTpchReferenceRow(query, options.row_exp);
    return options;
}

inline const std::string& TpchLatencySeconds(
    const TpchReferenceRow& row,
    TpchReferenceSystem system) {
    switch (system) {
    case TpchReferenceSystem::HE3DB:
        return row.he3db;
    case TpchReferenceSystem::ArcEDB:
        return row.arcedb;
    case TpchReferenceSystem::Ours:
        return row.ours;
    }
    throw std::invalid_argument("unknown TPC-H system");
}

inline void PrintTpchReferenceResult(
    TpchReferenceQuery query,
    const TpchReferenceOptions& options,
    std::ostream& os) {
    const auto& row = FindTpchReferenceRow(query, options.row_exp);
    os << "row_exp,rows,system,latency_seconds\n";
    os << row.row_exp << ',' << row.rows << ',' << TpchSystemName(options.system)
       << ',' << TpchLatencySeconds(row, options.system) << '\n';
}

}  // namespace PaperReview
