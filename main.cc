#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

#define int long long

const int INF = 1e18;

struct POINT{
    int pos_x, pos_y;
    string name;
};

struct CLUSTER{
    vector<pair<POINT,POINT>> path;
};

pair<vector<POINT>,vector<POINT>> in_inpout() {
    ifstream fin ("input.def", ifstream::in);
    fin.get();
    string s;
    bool start = false;
    while (!start && getline(fin, s)) {
        //cout << s << ' ';
        int other = 0, p_c = 0;
        for (int j = 0; j < s.size(); ++j) {
            //cout << j << ' ';
            if (s[j] == ';') ++p_c;
            if (s[j] != ' ') ++other;
        }
        if (other == 1 and p_c == 1)
            start = true;
    }
    vector<POINT> drivers;
    for (int i = 0; i < 32; ++i) {
        string name, b;
        fin >> b >> name;
        for (int j = 0; j < 23; ++j) fin >> b;
        int x, y;
        fin >> x >> y >> b >> b >> b;
        drivers.push_back({x, y, name});
    }
    string name, b;
    int x, y;
    vector<POINT> chip;
    while (fin >> name) {
        fin >> b >> b >> b >> b >> x >> y >> b >> b;
        chip.push_back({x, y, name});
    }
    return {drivers, chip};
}

int distance(POINT a, POINT b){
	return abs(a.pos_x-b.pos_x) + abs(a.pos_y-b.pos_y);
}

vector<POINT> getCluster(vector<POINT> p, int precision){
	int n = p.size();
	vector<POINT> center(16);
	for (int i = 0; i < 16 and i < p.size(); i++) center[i] = p[i];

	for (int b = 0; b < precision; b++) {
		vector<int> id(n);
		for (int i = 0; i < n; i++) {
			int MIN = 0, val = distance(p[i], center[0]);
			for (int j = 1; j < 16; j++) {
				if (distance(p[i], center[j]) <  val) {
					MIN = j;
					val = distance(p[i], center[j]);
				}
			}
			id[i] = MIN;
		}

		vector<vector<int>> X(16), Y(16); //coordinates x,y per cluter
		for (int i = 0; i < n; i++) {
			X[ id[i] ].push_back( p[i].pos_x );
			Y[ id[i] ].push_back( p[i].pos_y );
		}
		//CHANGE O(nlogn) --> O(n) (Amortized)
		for (int i = 0; i < 16; i++) nth_element(X[i].begin(), X[i].begin() + X[i].size()/ 2, X[i].end());
		for (int i = 0; i < 16; i++) nth_element(X[i].begin(), X[i].begin() + (X[i].size()+1)/ 2, X[i].end());
		for (int i = 0; i < 16; i++) nth_element(Y[i].begin(), Y[i].begin() + Y[i].size()/ 2, Y[i].end());
		for (int i = 0; i < 16; i++) nth_element(Y[i].begin(), Y[i].begin() + (Y[i].size()+1)/ 2, Y[i].end());
		//~ for (int i = 0; i < 16; i++) sort(X[i].begin(), X[i].end()), sort(Y[i].begin(), Y[i].end());
		//CHANGE
		vector<POINT> C(16); //center' of cluter using absolute difference formula
		for (int i = 0; i < 16; i++){
			C[i].pos_x = X[i][ X[i].size()/2 ];
			C[i].pos_y = Y[i][ Y[i].size()/2 ];
		}

		vector<int> DIST(16, INF); //This track the current nearest point to center'
		for (int i = 0; i < n; i++){
			if ( distance(p[i], C[ id[i] ]) < DIST[ id[i] ] ) {
				DIST[ id[i] ] = distance(p[i], C[ id[i] ]);
				center[ id[i] ] = p[i];
			}
		}
	}
	return center;
}

bool cmp_normal(pair<POINT, int>& a, pair<POINT, int>& b) {
    if (a.first.pos_x != b.first.pos_x) return a.first.pos_x < b.first.pos_x;
    return a.first.pos_y < b.first.pos_y;
}

vector<POINT> assign_inoutdrivers(vector<POINT>& center, vector<POINT>& in_out) {
    vector<POINT> assign_inout(32);
    //for (int i = 0; i < 32; ++i) assign_inout[i] = in_out[i];
    vector<pair<POINT, int>> ordenacion_centros;
    vector<pair<POINT, int>> ordenacion_in;
    vector<pair<POINT, int>> ordenacion_out;
    for (int i = 0; i < center.size(); ++i)
        ordenacion_centros.push_back({center[i], i});
    for (int i = 0; i < 16; ++i)
        ordenacion_in.push_back({in_out[i], i});
    for (int i = 0; i < 16; ++i)
        ordenacion_out.push_back({in_out[i+16], i+16});

    sort(ordenacion_centros.begin(), ordenacion_centros.end(), cmp_normal);
    sort(ordenacion_in.begin(), ordenacion_in.end(), cmp_normal);
    sort(ordenacion_out.begin(), ordenacion_out.end(), cmp_normal);

    for (int i = 0; i < 16; ++i) {
        assign_inout[ordenacion_centros[i].second] = ordenacion_in[i].first;
        assign_inout[ordenacion_centros[i].second + 16] = ordenacion_out[i].first;
    }
    return assign_inout;
}

pair<vector<vector<POINT>>, vector<POINT>> cluster_make(vector<POINT>& p, vector<POINT>& in_out) {
    int n = p.size();
    int precision = 20; // more precision, slower, more precise
	vector<POINT> center = getCluster(p, precision);

	vector<int> id(n);
	for (int i = 0; i < n; i++) {
		int MIN = 0, val = distance(p[i], center[0]);
		for (int j = 1; j < 16; j++) {
			if (distance(p[i], center[j]) <  val) {
				MIN = j;
				val = distance(p[i], center[j]);
			}
		}
		id[i] = MIN;
	}
	vector<vector<POINT>> sol(16, vector<POINT>());
	for (int i = 0; i < n; ++i) {
        sol[id[i]].push_back(p[i]);
	}

	return {sol, assign_inoutdrivers(center, in_out)};
}

int sz;

bool cmp (POINT& a, POINT& b) {
    if (a.pos_y/sz != b.pos_y/sz) return a.pos_y < b.pos_y;
    return a.pos_x < b.pos_x;
}

bool cmp2 (POINT& a, POINT& b) {
    if (a.pos_x/sz != b.pos_x/sz) return a.pos_x < b.pos_x;
    return a.pos_y < b.pos_y;
}

///Crear posibles funciones para ordenar

pair<int,CLUSTER> Mo_hor(vector<POINT>& l, int it, vector<POINT>& assignin_out, int sz_x) {
    sz = sz_x;
    int n = l.size();
    pair<int, CLUSTER> ans;
    sort(l.begin(), l.end(), cmp);
    CLUSTER sol;
    sol.path.push_back({assignin_out[it], l[0]});
    //cout << assignin_out[it].pos_y << ' ';
    long long man_dist = 0;
    for (int i = 1; i < n; ++i) {
        sol.path.push_back({l[i-1], l[i]});
        man_dist += abs (l[i].pos_x - l[i-1].pos_x) + abs(l[i].pos_y - l[i-1].pos_y);
    }
    sol.path.push_back({l[n-1], assignin_out[it+16]});
    //cout << assignin_out[it+16].pos_y << ' ';
    return {man_dist, sol};
}

pair<int,CLUSTER> Mo_ver(vector<POINT>& l, int it, vector<POINT>& assignin_out, int sz_y) {
    sz = sz_y;
    int n = l.size();
    pair<int, CLUSTER> ans;
    sort(l.begin(), l.end(), cmp2);
    CLUSTER sol;
    sol.path.push_back({assignin_out[it], l[0]});
    long long man_dist = 0;
    for (int i = 1; i < n; ++i) {
        sol.path.push_back({l[i-1], l[i]});
        man_dist += abs (l[i].pos_x - l[i-1].pos_x) + abs(l[i].pos_y - l[i-1].pos_y);
    }
    //cout << l[n-1].pos_x << ' ' << l[n-1].pos_y << '\n';
    sol.path.push_back({l[n-1], assignin_out[it+16]});
    return {man_dist, sol};
}

/**
Crear funciones para comparar valores
**/



pair<int,vector<CLUSTER>> best(vector<vector<POINT>>& l, vector<POINT>& assignin_out) {
    pair<int, vector<CLUSTER>> ans = {0, {}};
    for (int i = 0; i < 16; ++i) {
        pair<int, CLUSTER> best_option = {1e18, {}};
        int sz_x, sz_y;

        int min_x = INF, max_x = -INF, min_y = INF, max_y = -INF;
        int n = l.size();
        for (int j = 0; j < l[i].size(); ++j) {
            min_x = min(min_x, l[i][j].pos_x);
            max_x = max(max_x, l[i][j].pos_x);
            min_y = min(min_y, l[i][j].pos_y);
            max_y = max(max_y, l[i][j].pos_y);
        }
        sz_x = sqrt(max_x - min_x + 1);
        sz_y = sqrt(max_y - min_y + 1);
        int x = l[i].size();
        long double log_n = 0.001 / (1024);
        while (x > 0) {
            log_n *= 2;
            x /= 2;
        }
        for (long double j = 1/1.1; j < 1.1; j += log_n) {
            pair<int,CLUSTER> brute = Mo_hor (l[i], i, assignin_out, sz_x*j);
            pair<int,CLUSTER> s = Mo_ver (l[i], i, assignin_out, sz_y*j);
            if (brute.first < best_option.first) {
                best_option = brute;
            }
            else if (s.first < best_option.first) {
                best_option = s;
            }
        }
        ans.first += best_option.first;
        ans.second.push_back(best_option.second);
    }
    return ans;
}

void print_OUTPUT(const vector<CLUSTER>& l) {
    ofstream fout ("out.txt", ofstream::out);
    int n = l.size();
    for (int i = 0; i < 16; ++i) {
        int m = l[i].path.size();
        // fout << m+1 << ' ';
        for (int j = 0; j < m; ++j) {
            fout << "- BOGUS NET NAME\n";
            fout << "  ( " << l[i].path[j].first.name << " conn_in )\n";
            fout << "  ( " << l[i].path[j].second.name << " conn_out )\n";
            fout << ";\n";
            // fout << l[i].path[j].first.pos_x << ' ' << l[i].path[j].first.pos_y << ' ';
            //cout << l[i].path[0].first.pos_x << ' ' << l[i].path[0].first.pos_y << '\n';
        }
        // fout << l[i].path[m-1].second.pos_x << ' ' << l[i].path[m-1].second.pos_y << '\n';
        // fout << '\n';
    }
}

void path_making(vector<vector<POINT>>& l, vector<POINT>& assignin_out) {
    pair<int,vector<CLUSTER>> ans = best(l, assignin_out);
    print_OUTPUT(ans.second);
}

signed main() {
    pair<vector<POINT>, vector<POINT>> entrada = in_inpout();
    vector<POINT> in_out = entrada.first;
    vector<POINT> chip = entrada.second;
    pair<vector<vector<POINT>>, vector<POINT>> cluster_inout = cluster_make(chip, in_out);
    vector<vector<POINT>> cluster = cluster_inout.first;
    vector<POINT> assign_inout = cluster_inout.second;
    path_making(cluster, assign_inout);
}
