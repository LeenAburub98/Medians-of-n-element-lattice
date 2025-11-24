#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_set>

using namespace std;

// Function to find the join (least upper bound) of two elements in the lattice
int findJoin(const vector<vector<int>>& Matrix, int a, int b) {
    for (int i = 0; i < Matrix.size(); ++i) {
        if (Matrix[a][i] && Matrix[b][i]) {
            return i;
        }
    }
    return -1;  // No join found (shouldn't happen in a lattice)
}

// Function to find the meet (greatest lower bound) of two elements in the lattice
int findMeet(const vector<vector<int>>& Matrix, int a, int b) {
    for (int i = Matrix.size() - 1; i >= 0; --i) {
        if (Matrix[i][a] && Matrix[i][b]) {
            return i;
        }
    }
    return -1;  // No meet found (shouldn't happen in a lattice)
}

// Function to compare two elements (triples) in TL with multiple conditions
bool comparablelessthan(const vector<vector<int>>& Matrix, const vector<int>& t1, const vector<int>& t2, string& reason) {
    if ((Matrix[t1[0]][t2[0]] == 1 && Matrix[t1[1]][t2[1]] == 1 && Matrix[t1[2]][t2[2]] == 1) ||
        (Matrix[t1[0]][t2[0]] == 1 && Matrix[t1[1]][t2[2]] == 1 && Matrix[t1[2]][t2[1]] == 1) ||
        (Matrix[t1[0]][t2[1]] == 1 && Matrix[t1[1]][t2[0]] == 1 && Matrix[t1[2]][t2[2]] == 1) ||
        (Matrix[t1[0]][t2[1]] == 1 && Matrix[t1[1]][t2[2]] == 1 && Matrix[t1[2]][t2[0]] == 1) ||
        (Matrix[t1[0]][t2[2]] == 1 && Matrix[t1[1]][t2[0]] == 1 && Matrix[t1[2]][t2[1]] == 1) ||
        (Matrix[t1[0]][t2[2]] == 1 && Matrix[t1[1]][t2[1]] == 1 && Matrix[t1[2]][t2[0]] == 1)) {
        reason = "true";
        return true;
    } else {
        reason = "false.";
        return false;
    }
}

int main() {
    int n;
    cout << "Enter the number of elements in the lattice: ";
    cin >> n;

    vector<vector<int>> Matrix(n, vector<int>(n));

    cout << "Enter the adjacency matrix of the lattice:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> Matrix[i][j];
        }
    }

    ofstream outFile("lattice_output.txt");

    // Writing meet and join for pairs
   // outFile << "Meet(i, j) | Join(i, j)\n";
    //outFile << "-------------------------------\n";
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int join = findJoin(Matrix, i, j);
            int meet = findMeet(Matrix, i, j);
            //outFile << "Meet(" << i << ", " << j << ") = " << meet
              //      << " | Join(" << i << ", " << j << ") = " << join << "\n";
        }
    }

   /* outFile << " x M ( y J ( x M z))";
    outFile << "-------------------------------\n";*/
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
                for(int k = j + 1; k < n; ++k){
            int joinij = findJoin(Matrix, i, j);
            int meetij = findMeet(Matrix, i, j);
            int joinik = findJoin(Matrix, i, k);
            int meetik = findMeet(Matrix, i, k);
            int joinjk = findJoin(Matrix, j, k);
            int meetjk = findMeet(Matrix, j, k);

            int jJMik = findJoin(Matrix, j, meetik);
            int iJMjk = findJoin(Matrix, i, meetjk);
            int kJMij = findJoin(Matrix, k, meetij);
            int iMjJMik = findMeet(Matrix, i, jJMik);
            int iMkJMij = findMeet(Matrix, i, kJMij);
            int jMiJMjk = findMeet(Matrix, j, iJMjk);
            int jMkJMji = findMeet(Matrix, j, kJMij);
            int kMiJMjk = findMeet(Matrix, k, iJMjk);
            int kMjJMik = findMeet(Matrix, k, jJMik);

            int totJoin = findJoin (Matrix, findJoin (Matrix, iMjJMik, iMkJMij), findJoin(Matrix, findJoin (Matrix, jMiJMjk, jMkJMji), findJoin (Matrix, kMiJMjk, kMjJMik)));

            int jMJik = findMeet(Matrix, j, joinik);
            int iMJjk = findMeet(Matrix, i, joinjk);
            int kMJij = findMeet(Matrix, k, joinij);
            int iJjMJik = findJoin(Matrix, i, jMJik);
            int iJkMJij = findJoin(Matrix, i, kMJij);
            int jJiMJjk = findJoin(Matrix, j, iMJjk);
            int jJkMJji = findJoin(Matrix, j, kMJij);
            int kJiMJjk = findJoin(Matrix, k, iMJjk);
            int kJjMJik = findJoin(Matrix, k, jMJik);

            int totMeet = findMeet (Matrix, findMeet (Matrix, iJjMJik, iJkMJij), findMeet(Matrix, findMeet (Matrix, jJiMJjk, jJkMJji), findMeet (Matrix, kJiMJjk, kJjMJik)));



            outFile << "studying ij ik jk: "<< "\n";
            outFile << "........................ " << i << j << k << "........................ \n";
            outFile << " meet " << i << j << " is " << meetij << "\n meet " << i << k << " is " << meetik << "\n meet " << j << k << " is " << meetjk << "\n join " << i << j << " is " << joinij
            << "\n join " << i << k << " is " << joinik << "\n join " << j << k << " is " << joinjk << "\n";

            outFile << i << " join (" << j << " meet " << k << ") is " << iJMjk << "\n";
            outFile << j << " join (" << i << " meet " << k << ") is " << jJMik << "\n";
            outFile << k << " join (" << i << " meet " << j << ") is " << kJMij << "\n";

            outFile << i << " meet (" << j << " join " << k << ") is " << iMJjk << "\n";
            outFile << j << " meet (" << i << " join " << k << ") is " << jMJik << "\n";
            outFile << k << " meet (" << i << " join " << j << ") is " << kMJij << "\n";


            outFile << i << " join (" << j << " meet (" << i << " join " << k << ")) is " << iJjMJik << "\n" ;
            outFile << i << " join (" << k << " meet (" << i << " join " << j << ")) is " << iJkMJij << "\n" ;

            outFile << j << " join (" << i << " meet (" << j << " join " << k << ")) is " << jJiMJjk << "\n" ;
            outFile << j << " join (" << k << " meet (" << j << " join " << i << ")) is " << jJkMJji << "\n" ;

            outFile << k << " join (" << i << " meet (" << k << " join " << j << ")) is " << kJiMJjk << "\n" ;
            outFile << k << " join (" << j << " meet (" << k << " join " << j << ")) is " << kJjMJik << "\n" ;

            outFile << "total meet of the above 6 terms: " << totMeet << "\n" ;

            outFile << i << " meet (" << j << " join (" << i << " meet " << k << ")) is " << iMjJMik << "\n" ;
            outFile << i << " meet (" << k << " join (" << i << " meet " << j << ")) is " << iMkJMij << "\n" ;

            outFile << j << " meet (" << i << " join (" << j << " meet " << k << ")) is " << jMiJMjk << "\n" ;
            outFile << j << " meet (" << k << " join (" << j << " meet " << i << ")) is " << jMkJMji << "\n" ;

            outFile << k << " meet (" << i << " join (" << k << " meet " << j << ")) is " << kMiJMjk << "\n" ;
            outFile << k << " meet (" << j << " join (" << k << " meet " << j << ")) is " << kMjJMik << "\n" ;

            outFile << "total join of the above 6 terms: " << totJoin << "\n" ;

                }
        }
    }

    // Automatically find TL set from the matrix based on the condition
    vector<vector<int>> TL;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                int meet_ij = findMeet(Matrix, i, j);
                int meet_ik = findMeet(Matrix, i, k);
                int meet_jk = findMeet(Matrix, j, k);

                int join_of_meets = findJoin(Matrix, meet_ij, findJoin(Matrix, meet_ik, meet_jk));

                int join_ij = findJoin(Matrix, i, j);
                int join_ik = findJoin(Matrix, i, k);
                int join_jk = findJoin(Matrix, j, k);

                int meet_of_joins = findMeet(Matrix, join_ij, findMeet(Matrix, join_ik, join_jk));

                if (join_of_meets != meet_of_joins) {
                    TL.push_back({i, j, k});
                }
            }
        }
    }

    int m = TL.size();

    // Vector storing the PI size for each element in TL
    vector<int> ps(m, 0);

    // Initialize comparison table and PI table
    vector<vector<int>> ct(m, vector<int>(m, 0));  // Comparison table
    vector<vector<int>> pe(m, vector<int>(n, 0));  // Permitted intervals table

    // Fill in the comparison table (ct)
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            string reason;
            if (comparablelessthan(Matrix, TL[i], TL[j], reason)) {
                ct[i][j] = 1;
            }
        }
    }

    // Fill in the permitted intervals for each element in TL
    for (int i = 0; i < m; ++i) {
        const auto& triple = TL[i];
        int a = triple[0], b = triple[1], c = triple[2];

        // Calculate meets and joins for the current triple
        int meet_ab = findMeet(Matrix, a, b);
        int meet_ac = findMeet(Matrix, a, c);
        int meet_bc = findMeet(Matrix, b, c);
        int join_of_meets = findJoin(Matrix, meet_ab, findJoin(Matrix, meet_ac, meet_bc));
        int join_ab = findJoin(Matrix, a, b);
        int join_ac = findJoin(Matrix, a, c);
        int join_bc = findJoin(Matrix, b, c);
        int meet_of_joins = findMeet(Matrix, join_ab, findMeet(Matrix, join_ac, join_bc));

        // Populate the PI set and its size
        for (int x = 0; x < n; ++x) {
            if (Matrix[join_of_meets][x] == 1 && Matrix[x][meet_of_joins] == 1) {
                pe[i][x] = 1;  // x is in the PI of vector[i]
                ps[i]++;
            }
        }
    }

    // Output the PI sizes
    outFile << "\nPermitted Intervals Sizes for TL elements:\n";
    for (int i = 0; i < m; ++i) {
        outFile << "Size of PI for element " << i << " (" << TL[i][0] << "," << TL[i][1] << "," << TL[i][2] << "): " << ps[i] << "\n";
    }

    // Output the permitted intervals
    outFile << "\nPermitted Intervals for TL elements:\n";
    for (int i = 0; i < m; ++i) {
        outFile << "PI(" << i << ") = { ";
        for (int j = 0; j < n; ++j) {
            if (pe[i][j] == 1) {
                outFile << j << " ";
            }
        }
        outFile << "}\n";
    }

    // Begin processing combinations
    long long s = 1;
    for (int t = 0; t < m; ++t) {
        s *= ps[t];  // Calculate s based on the size of PIs
    }

    vector<bool> truemedian(s, true);
    vector<vector<int>> med(s, vector<int>(m));

    for (long long c = 0; c < s; ++c) {
        long long temp = c;
        vector<int> h(m);
        for (int t = 0; t < m; ++t) {
            h[t] = temp % ps[t];
            temp /= ps[t];
        }

        for (int t = 0; t < m; ++t) {
            int d = h[t];
            for (int j = 0; j < n; ++j) {
                if (pe[t][j] == 1) {
                    if (d == 0) {
                        med[c][t] = j;
                        break;
                    }
                    d--;
                }
            }
        }

        // Validate the truemedian
        for (int ia = 0; ia < m; ++ia) {
            for (int ib = 0; ib < m; ++ib) {
                if (ct[ia][ib] == 1) {
                    if (Matrix[med[c][ia]][med[c][ib]] == 0) {
                        truemedian[c] = false;
                        break;
                    }
                }
            }
        }
    }

    int om = 0;
    for (long long c = 0; c < s; ++c) {
        if (truemedian[c]) {
            om++;
        }
    }

    outFile << "\nTotal number of true medians: " << om << "\n";

    outFile.close();
    cout << "Output written to lattice_output.txt\n";

    return 0;
}
