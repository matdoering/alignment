#include "Alignment.h"
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>

namespace Align{


EditDistance::EditDistance(NucleicAcidType nucleicAcidType) {
    switch(nucleicAcidType) {
        case NucleicAcidType::DNA:
            initDNADistances();
            break;
        case NucleicAcidType::RNA:
            initRNADistances();
            break;
        default:
            abort(); // Unknown nucleic acid type
    }
}

void EditDistance::initRNADistances() {
    m_distances["--"] = 0.;
    m_distances["-A"] = -2.;
    m_distances["-C"] = -2.;
    m_distances["-G"] = -2.;
    m_distances["-U"] = -2.;

    m_distances["AA"] = 3.0;
    m_distances["AC"] =  -1.;
    m_distances["AG"] =  -1.;
    m_distances["AU"] = -1.;
    m_distances["A-"] = -2.;

    m_distances["CA"] = -1.;
    m_distances["CC"] = 3.;
    m_distances["CG"] = -1.;
    m_distances["CU"] = -1.;
    m_distances["C-"] = -2.;

    m_distances["GA"] = -1.;
    m_distances["GC"] = -1.;
    m_distances["GG"] = 3.;
    m_distances["GU"] = -1.;
    m_distances["G-"] = -2.;

    m_distances["UA"] =  -1.;
    m_distances["UC"] = -1.;
    m_distances["UG"] = -1.;
    m_distances["UU"] =  3.;
    m_distances["U-"] = -2.;
}

void EditDistance::initDNADistances() {
    m_distances["--"] = 0.;

    m_distances["-A"] = -2.;
    m_distances["-C"] = -2.;
    m_distances["-G"] = -2.;
    m_distances["-T"] = -2.;

    m_distances["AA"] = 3.0;
    m_distances["AC"] =  -1.;
    m_distances["AG"] =  -1.;
    m_distances["AT"] = -1.;
    m_distances["A-"] = -2.;

    m_distances["CA"] = -1.;
    m_distances["CC"] = 3.;
    m_distances["CG"] = -1.;
    m_distances["CT"] = -1.;
    m_distances["C-"] = -2.;

    m_distances["GA"] = -1.;
    m_distances["GC"] = -1.;
    m_distances["GG"] = 3.;
    m_distances["GT"] = -1.;
    m_distances["G-"] = -2.;

    m_distances["TA"] =  -1.;
    m_distances["TC"] = -1.;
    m_distances["TG"] = -1.;
    m_distances["TT"] =  3.;
    m_distances["T-"] = -2.;
}

float EditDistance::dist(const char& c1, const char& c2) const {
    std::string key = std::string(1, c1) + std::string(1, c2);
    return m_distances.at(key);
}

namespace{
    size_t getMaxIdx(std::vector<float> vals) {
        size_t maxIdx = 0;
        for (size_t i =0; i < vals.size(); ++i) {
            if (vals[i] > vals[maxIdx]) { 
                maxIdx = i;
            }
        }
        return maxIdx;
    }

    void setMaxValue(matrix& M, int row, int col, const EditDistance& dist, const char& s1, const char& s2) {
        /// Sets the maximum value in the matrix M
        float leftValue = M[row][col-1].value + dist.dist(s_gapChar, s2);
        float topValue = M[row-1][col].value + dist.dist(s1, s_gapChar);
        float diagValue = M[row-1][col-1].value + dist.dist(s1, s2);
        std::vector<float> data = {leftValue, topValue, diagValue};

        // determine max
        auto maxIdx = getMaxIdx(data);
        Elem elem;
        BacktrackingChoice btChoice = BacktrackingChoice::INVALID;
        switch(maxIdx) {
            case 0:
                btChoice = BacktrackingChoice::LEFT;
                elem = Elem(data[maxIdx], std::make_pair(row, col-1), btChoice);
                break;
            case 1:
                btChoice = BacktrackingChoice::TOP;
                elem = Elem(data[maxIdx], std::make_pair(row-1, col), btChoice);
                break;
            case 2:
                btChoice = BacktrackingChoice::DIAG;
                elem = Elem(data[maxIdx], std::make_pair(row-1, col-1), btChoice);
                break;
            default:
                abort(); // impossible max idx
        }
        M[row][col] = elem;
    }

}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (auto& row : m.M) {
        for (auto& elem : row) {
            os << elem.value << " ";
        }
        os << "\n";
    }
    return os;    
}

NucleicAcidType detectNucleicAcidType(const std::string& s) {
    size_t countT = 0;
    size_t countU = 0;
    for (auto& c : s) {
        if (c == 'T') {
            ++countT;
        } else if (c == 'U') {
            ++countU;
        }
    }
    if (countT > 0 && countU > 0) {
        throw std::invalid_argument("Input sequence used both T's and U's");
    }
    if (countU > 0) {
        // assume RNA
        return NucleicAcidType::RNA;
    } else {
        return NucleicAcidType::DNA;
    }
}

std::ostream& operator<<(std::ostream& os, const Alignment& ali) {
    os << "s1: " << ali.m_s1 << "\n"
       << "s2: " << ali.m_s2 << "\n"
       << "Alignment:" << "\n"
       << ali.m_s1A << "\n" 
       << ali.m_s2A << "\n"
       << "Score: " << ali.m_aliScore;
    return os;
}


void Alignment::align() {
    // detect nucleic acid type and choose edit distance
    NucleicAcidType s1Type = detectNucleicAcidType(m_s1);
    NucleicAcidType s2Type = detectNucleicAcidType(m_s2);
    if (s1Type != s2Type) {
        throw std::invalid_argument("Input sequence types (DNA/RNA) do not match");
    }
    EditDistance dist(s1Type);
    
    Matrix m(m_s1.size()+1, m_s2.size()+1);
    initMatrix(m, dist);
    fillMatrix(m, dist);
    traceback(m);
}

void Alignment::fillMatrix(Matrix&m, const EditDistance& dist) {
    // fill matrix
    for (size_t row = 1; row < m.nrow; ++row) {
        for (size_t col = 1; col < m.ncol; ++col) {
            setMaxValue(m.M, row, col, dist, m_s1[row-1], m_s2[col-1]);
        }
    }
}

void Alignment::initMatrix(Matrix& m, const EditDistance& dist) {
    matrix& M = m.M;
    BacktrackingChoice invalidChoice = BacktrackingChoice::DIAG;

    // set initial element
    M[0][0] = Elem(0, std::make_pair(-1, -1), invalidChoice);

    // init column 0
    for (size_t row = 1; row < m.nrow; ++row) {
        float val = M[row-1][0].value + dist.dist(m_s1[row-1], s_gapChar);
        M[row][0] = Elem(val, std::make_pair(-1, row-1), invalidChoice);
    }
    // init row 0
    for (size_t col = 1; col < m.ncol; ++col) {
        float val = M[0][col-1].value + dist.dist(s_gapChar, m_s2[col-1]);
        M[0][col] = Elem(val, std::make_pair(col-1, -1), invalidChoice);
    }
}

void Alignment::traceback(const Matrix& m) {
    const matrix& M = m.M;
    std::pair<int,int> btCoord = std::make_pair(m.nrow-1,m.ncol-1);
    std::pair<int,int> stopCoord = std::make_pair(-1,-1);

    size_t expectedAliSize = m.nrow > m.ncol ? m.nrow-1 : m.ncol-1;
    m_s1A = "";
    m_s1A.resize(expectedAliSize);
    m_s2A = "";
    m_s2A.resize(expectedAliSize);

    m_aliScore = 0.0;
    while (btCoord != stopCoord) {
        const Elem& elem = M[btCoord.first][btCoord.second];
        // build alignment string
        BacktrackingChoice btChoice = elem.backtrackingChoice;
        const char& c1 = m_s1[btCoord.first];
        const char& c2 = m_s2[btCoord.second];
        switch(btChoice) {
            case BacktrackingChoice::DIAG:
                m_s1A += c1;
                m_s2A += c2;
                break;
            case BacktrackingChoice::LEFT:
                m_s1A += "-";
                m_s2A += c2;
                break;
            case BacktrackingChoice::TOP:
                m_s1A += c1;
                m_s2A += "-";
                break;
            case BacktrackingChoice::INVALID:
                abort();
        }
        m_aliScore += M[btCoord.first][btCoord.second].value;
        btCoord = elem.backtrackingCoordinate;
    }
    std::reverse(m_s1A.begin(), m_s1A.end());
    std::reverse(m_s2A.begin(), m_s2A.end());
}

}


