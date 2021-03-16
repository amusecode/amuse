/*MIT License

Copyright (c) 2021 Keigo Nitadori

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

/*
-added support for QD in addition to DD by Thomas Schano  QD is available at http://crd-legacy.lbl.gov/~dhbailey/mpdist/ (https://github.com/GFTwrt/QD)
*/

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
