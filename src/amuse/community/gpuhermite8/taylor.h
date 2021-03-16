template<
	int      beg,
	int      end,
	typename real_t,
	typename vect_t>
struct TaylorInner{
	static vect_t eval(
	   const real_t h,
	   const vect_t y[])
	{
		return y[beg] + (h * (real_t(1)/real_t(beg+1)))
			* TaylorInner<beg+1, end, real_t, vect_t>::eval(h, y);
			
	}
};

template<
	int      end,
	typename real_t,
	typename vect_t>
struct TaylorInner <end, end, real_t, vect_t>{
	static vect_t eval(
	   const real_t h,
	   const vect_t y[])
	{
		return y[end];
			
	}
};

template<
	int      p,
	typename real_t,
	typename vect_t>
static vect_t taylor(
	   const real_t h,
	   const vect_t y[])
{
	return TaylorInner<0, p, real_t, vect_t>::eval(h, y);
}
