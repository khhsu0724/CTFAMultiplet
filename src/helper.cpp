#include "helper.hpp"

size_t ed::choose(size_t n, size_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;  // Use symmetry

    size_t res = 1;
    for (size_t i = 1; i <= k; ++i) {
        res = res * (n - k + i) / i;
    }
    return res;
}

bool ed::is_pw2(ulli x) {return !(x == 0) && !(x & (x - 1));}

ulli ed::next_perm(ulli v) {
	ulli t = v | (v - 1);
	return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(v) + 1));
}

// Enumarate states that include bits for "inc" variable using recursion
void ed::enum_states(std::vector<ulli>& states, ulli n, ulli k, ulli inc, ulli s) {
	if (k > n) return;
	if (k == 0 && n < (__builtin_ctzll(inc)+1)) {
		states.emplace_back(s);
		return;
	}
	enum_states(states,n-1,k-1,inc,s+(1<<(n-1))); // bit is 1
	if (!(inc & (1<<(n-1)))) enum_states(states,n-1,k,inc,s); // bit is 0
}

// Add valence state and core state 
ulli ed::add_bits(ulli b1, ulli b2, int b1size, int b2size) {
	ulli b1d = b1 & ((BIG1 << (b1size/2)) - 1);
	ulli b2d = b2 & ((BIG1 << (b2size/2)) - 1);
	return (b1^b1d) << b2size | (b2^b2d) << (b1size/2) | b1d << (b2size/2) | b2d;
}

int ed::count_bits(ulli b) {
	int c;
	for (c = 0; b; c++) b &= b - 1;
	return c;
}

vecc ed::vec_conj(vecc vin) {
	vecc vout(vin.size(),0);
	#pragma omp parallel for
	for (int i = 0; i < vin.size(); ++i) vout[i] = std::conj(vin[i]);
	return vout;
}

vecc ed::ctranspose(const vecc& mat, size_t m, size_t n) {
	// Conjugate Transpose of a matrix
	vecc trans(m*n,0);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < m; j++) {
			trans.at(j*n+i) = std::conj(mat.at(i*m+j));
		}
	}
	return trans;
}


std::vector<int> ed::distribute(int num_h, int num_at) {
	int base = num_h / num_at;
	int mod = num_h % num_at;
	std::vector<int> dist(num_at,base);
	for (int i = 0; i < mod; ++i) ++dist[i];
	return dist;
}

std::string ed::format_duration(std::chrono::milliseconds ms) {
    using namespace std::chrono;
    auto secs = duration_cast<seconds>(ms);
    ms -= duration_cast<milliseconds>(secs);
    auto mins = duration_cast<minutes>(secs);
    secs -= duration_cast<seconds>(mins);
    auto hour = duration_cast<hours>(mins);
    mins -= duration_cast<minutes>(hour);

    std::stringstream ss;
    ss << hour.count() << " Hours : " << mins.count() << " Minutes : " << secs.count() << " Seconds : " << ms.count() << " Milliseconds";
    return ss.str();
}

void ed::print_progress(double frac, double all) {
	// Print at least every 2%
	if (((int)frac+1)%((int)all/50+1) && frac != all) return;
	double progress = frac/all;
    int barWidth = 40;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
	if (progress == 1) std::cout << std::endl;
	return;
}

void ed::parse_num(std::string complex_string, dcomp& complex_num) {
	// Parse complex number from string. Format: (x,y)
	complex_string = complex_string.substr(complex_string.find("(")+1);
	complex_string = complex_string.substr(0,complex_string.find(")"));
	size_t comma =  complex_string.find(",");
	double real = stod(complex_string.substr(0,comma));
	double imag = stod(complex_string.substr(comma+1));
	complex_num = dcomp(real,imag);
	return;
}

void ed::write_vecc(vecc vec, size_t x, size_t y, std::string file_dir, std::string delim) {
	// Specifically write to python readable format
	std::ofstream matfile;
    matfile.open(file_dir);
	for (int j = 0; j < y; ++j) {
    	for (int i = 0; i < x; ++i) {
    		double r = std::real(vec.at(i+x*j));
    		double c = std::imag(vec.at(i+x*j));
    		matfile << std::setprecision(5) << r;
    		if (c >= 0) matfile << "+";
    		matfile << std::setprecision(5) << c << "j";
    		if (i < x-1) matfile << delim;
    	}
    	matfile << "\n";
    }
    matfile.close();
};
