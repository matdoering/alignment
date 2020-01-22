#pragma once

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <utility>

namespace Align {

    static const char s_gapChar = '-';

enum class NucleicAcidType {
        DNA,
        RNA,
        DNA_OR_RNA
    };


    struct EditDistance{
        public:
            EditDistance(NucleicAcidType naType);
            float dist(const char& c1, const char& c2) const;
        private:
        void initDNADistances();
        void initRNADistances();
        std::map<std::string, float> m_distances;
            /// Edit distances
    };
    
    enum class BacktrackingChoice {
        LEFT = 1, // move from left cell in the matrix
        TOP, // move from top cell in the matrix
        DIAG, // move from diagonally top left cell in the matrix
        INVALID
    };

    struct Elem{
        Elem() = default;
        Elem(float val, std::pair<int,int>&& btCoord, BacktrackingChoice btChoice) : value(val), backtrackingCoordinate(btCoord), backtrackingChoice(btChoice)
        {}

        float value;
            /// Similarity value
        std::pair<int,int> backtrackingCoordinate;
            /// The coordinate in the matrix from which we reached this Elem during backtracking
        BacktrackingChoice backtrackingChoice;
            /// The backtracking cell direction where we come from (diag / top / left)

    };

    using matrix = std::vector<std::vector<Elem>>;

    /*
     * Container for the dynamic programming matrix
     */
    struct Matrix {
        Matrix(size_t nrows, size_t ncols) : nrow(nrows), ncol(ncols), M(matrix(nrows, std::vector<Elem>(ncols))) {
        }
        matrix M;
        size_t nrow;
        size_t ncol;
        friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    };

    /*
     * Object representing an alignment
     */
    class Alignment{

        public:
        Alignment(std::string s1, std::string s2) : m_s1(s1), m_s2(s2) {
        }
        void align(bool isLocalAlignment);
        friend std::ostream& operator<<(std::ostream& os, const Alignment& ali);

        private:
            void traceback(const Matrix& M, bool isLocalAlignment);
            void initMatrixGlobal(Matrix& M, const EditDistance& dist);
                /// Initialize similarity matrix for global alignment
            void initMatrixLocal(Matrix& M, const EditDistance& dist);
                /// Initialize similarity matrix for local alignment
            void fillMatrix(Matrix& M, const EditDistance& dist);
                /// Fill the similarity matrix using dynamic programming
            std::string m_s1;
                // Input sequence s1
            std::string m_s2;
                // Input sequence s2
            std::string m_s1A;
                /// Aligned sequence s1
            std::string m_s2A;
                /// Aligned sequence s2
            float m_aliScore;
                /// Alignment score
                /// gap character used
    };
   
        
}
