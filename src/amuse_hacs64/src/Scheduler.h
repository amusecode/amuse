#ifndef _Scheduler_H_
#define _Scheduler_H_

#include <algorithm>
enum 
{
	RUNGMAX = 32, 
	REGFLAG = (1 << 31)
};


// IEEE complient ONLY !!!!!!!!!
static const int exp_of(const double dt) 
{
	union {double d; unsigned long long u;} pack;
	pack.d = dt;
	return (pack.u >> 52) & 0xfff;
};

struct Scheduler 
{
	double dtmax;
	double dt_tick;
	unsigned long long tsysU;
	int min_rung, max_rung;      // the timestep is dtmax/(1 << rung)
	std::vector<int> list[RUNGMAX];
	Scheduler() : dtmax(-1.0){}
	~Scheduler() {};

	Scheduler(const double _dtmax) : dtmax(_dtmax) {
		for (int i = 0; i < RUNGMAX; i++) {
			list[i].reserve(64);
			list[i].clear();
		}
		dt_tick = dtmax / (1LLU << (RUNGMAX - 1));
		fprintf(stderr, " dt_tick= %lg\n", dt_tick);
		tsysU = 0;
		min_rung = 0;
		max_rung = RUNGMAX - 1;
	}

  void flush() {flush_list();}
	void flush_list() {
		for (int i = 0; i < RUNGMAX; i++)
			list[i].clear();
	}

	const int move(const int rung_old, const int rung_new, const int idx) {
//    assert(0);
		std::vector<int>::iterator it = std::find(list[rung_old].begin(), list[rung_old].end(), idx);
		assert(it != list[rung_old].end());
		list[rung_old].erase(it);

		list[rung_new].push_back(idx);
		return rung_new;
	}

	void remove(const int rung, const int idx) {
		std::vector<int>::iterator it = std::find(list[rung].begin(), list[rung].end(), idx);
		assert(it != list[rung].end());
		list[rung].erase(it);
	}

	const int push_particle(const int index, const int rung0, const bool isreg = false) {
		const int rung = std::max(rung0, min_rung);
		assert(rung < RUNGMAX);
		list[rung].push_back(index | (isreg ? REGFLAG : 0));
		return rung;
	}

	const int push_particle(const int index, const double dt, const bool isreg = false) {
		const int rung = std::max(exp_of(dtmax) - exp_of(dt), min_rung);
		assert(rung < RUNGMAX);
		list[rung].push_back(index | (isreg ? REGFLAG : 0));
		return rung; //dtmax/(1U << rung);
	};

	const int get_rung(const double dt) {
		const int rung = std::max(exp_of(dtmax) - exp_of(dt), min_rung);
		assert(rung < RUNGMAX);
		return rung; 
	}
	const int get_rung(const int rung0)
  {
		const int rung = std::max(rung0, min_rung);
		assert(rung < RUNGMAX);
		return rung; 
	}

	const double get_dt(const int rung) {
		unsigned long long dtU = 1LLU << (RUNGMAX - 1 - rung);
		return dt_tick * dtU;
	}

	const double pull_active_list(std::vector<int> &irr_list, std::vector<int> &reg_list) {
		irr_list.clear();
		reg_list.clear();
		max_rung = 0;
		for (int i = 0; i < RUNGMAX; i++) {
			if (list[i].size() > 0) max_rung = i;
		}
		const unsigned long long dtU    = 1LLU << (RUNGMAX - 1 - max_rung);
    const unsigned long long tnextU = (tsysU/dtU + 1)*dtU;
		const unsigned long long flags  = tsysU ^ tnextU;

		for (int i = 0; i < RUNGMAX; i++) {
			if (flags & (1LLU << i)) {
				if (list[RUNGMAX - 1 - i].size() > 0) { 
					for (size_t j = 0; j < (const size_t)list[RUNGMAX - 1 - i].size(); j++) {
						const int idx = list[RUNGMAX - 1 - i][j];
						if ((idx & REGFLAG) == 0) irr_list.push_back(idx & ~REGFLAG);
						else                      reg_list.push_back(idx & ~REGFLAG);
					}
					list[RUNGMAX - 1 - i].clear();
				}
				min_rung = RUNGMAX - 1 - i;
			}
		}
		tsysU = tnextU;
		return dtmax / (1LLU << max_rung);
	}

	void debug_dump() const { 

		FILE *ferr = stderr;
		int nirr_tot = 0;
		int nreg_tot = 0;
		for (int i = 0; i < RUNGMAX; i++) {
			const int nirr = get_nirr(i);
			const int nreg = get_nreg(i);
			const int count = list[i].size();
			nirr_tot += nirr;
			nreg_tot += nreg;
			assert(nirr + nreg == count);
			fprintf(ferr, " %4d : %6d [ nirr= %6d   nreg= %6d ]  %g\n", i, count, nirr, nreg, dtmax/(1LLU<<i));
		}
		fprintf(ferr, " total: %6d [ nirr= %6d   nreg= %6d ]\n",  nirr_tot + nreg_tot, nirr_tot, nreg_tot);
	}

	const int get_nreg() const {
		int nreg = 0;
		for (int i = 0; i < RUNGMAX; i++)
			nreg += get_nreg(i);
		return nreg;
	}

	const int get_nirr() const {
		int nirr = 0;
		for (int i = 0; i < RUNGMAX; i++)
			nirr += get_nirr(i);
		return nirr;
	}
	const int get_nreg(const int i) const {
		int nreg = 0;
		const int count = list[i].size();
		for (int j = 0; j < count; j++) 
			nreg += (list[i][j] & REGFLAG) != 0;
		return nreg;
	}

	const int get_nirr(const int i) const {
		int nirr = 0;
		const int count = list[i].size();
		for (int j = 0; j < count; j++) 
			nirr += (list[i][j] & REGFLAG) == 0;
		return nirr;
	}

	const double get_tsys() const {
		return dt_tick * tsysU;
	}

	const double get_dt_pred(const int rung) const {
		const unsigned long long dtU = 1LLU << (RUNGMAX - 1 - rung);
		const unsigned long long tmp = (rung >= min_rung) ? dtU : tsysU & (dtU - 1);
		return dt_tick * tmp;
	}

	const double get_dt_corr(const int rung) const {
		return dtmax / (1LLU << rung) ;
	}
};

#endif // _Scheduler_H_ 
