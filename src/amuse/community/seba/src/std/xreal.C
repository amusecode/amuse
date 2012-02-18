
#include "stdinc.h"

#ifndef TOOLBOX

#ifdef USE_XREAL

inline int xadd(xfrac_t &a, xfrac_t b)
{
    a += b;
    return (a < b ? 1 : 0);	// a < b ==> overflow
}

inline int xsub(xfrac_t &a, xfrac_t b)
{
    xfrac_t aa = a;
    a -= b;
    return (a > aa ? 1 : 0);	// a > aa ==> underflow
}

// Constructors:

xreal::xreal() {i = f = 0;}					// default

//xreal::xreal(xreal x) {					// x = xreal
//    i = x.i;
//    f = x.f;
//}

xreal::xreal(xint_t ii, xfrac_t ff) {
    i = ii;
    f = ff;
}

xreal::xreal(int ii) {						// x = int
    i = ii;
    f = 0;
}

xreal::xreal(real x)						// x = real
{
    real ii = floor(x);			// split integer/fraction
    real ff = x - ii;

    // Check for over/underflow.

    if (ii <= -TWO63) {			// comparison with TWO63N
					// here fails!

	// i = TWO63N;			// does not work on Linux:
					// RHS is interpreted as int

	// * Ugly! *

	i = 1;
	i = i<<63;			// should work (overflow/wrap)...

	f = 0;				// smallest possible xreal

    } else if (ii >= TWO63) {

	// i = TWO63M;			// does not work on linux:
					// RHS is interpreted as int

	// * Ugly! *

	i = 1;
	i = i <<63;
	i--;				// should work (overflow/wrap)...

	f = 0; f--;			// largest possible xreal

    } else {
	i = (xint_t)ii;
	f = (xfrac_t)(TWO64*ff);
    }
}

real xreal::to_real()
{
    // Explicitly avoid obvious roundoff problems near zero.

    if (i == -1) {

	return -((0-f)*TWO64I);	// TWO64 <--> 0!
		
    } else
	return i + f*TWO64I;	// should be OK for other real values
}

// Unary -, binary +, -, +=, -= (don't check for over/underflow):

xreal xreal::operator - () {return xreal(-i-1, 0-f);}

xreal xreal::operator + (const xreal y) {
    xfrac_t sumf = f;
    xint_t sumi = i + y.i + xadd(sumf, y.f);
    return xreal(sumi, sumf);
}

xreal xreal::operator - (const xreal y) {
    xfrac_t sumf = f;
    xint_t sumi = i - y.i - xsub(sumf, y.f);
    return xreal(sumi, sumf);
}

xreal& xreal::operator += (const xreal y) {
    i += y.i + xadd(f, y.f);
    return *this;
}

xreal& xreal::operator -= (const xreal y) {
    i -= y.i + xsub(f, y.f);
    return *this;
}
xreal xreal::operator + (const real y) {
    return *this + (xreal)y;
}

xreal xreal::operator - (const real y) {
    return *this - (xreal)y;
}

xreal& xreal::operator += (const real y) {
    *this += (xreal)y;
    return *this;
}

xreal& xreal::operator -= (const real y) {
    *this -= (xreal)y;
    return *this;
}

// Logical xreal operators ==, !=, <, <=, >, >=:

bool xreal::operator == (const xreal y) {
    return (i == y.i && f == y.f);
}

bool xreal::operator != (const xreal y) {
    return (i != y.i || f != y.f);
}

bool xreal::operator < (const xreal y) {
    return (i < y.i || (i == y.i && f < y.f));
}

bool xreal::operator <= (const xreal y) {
    return (i < y.i || (i == y.i && f <= y.f));
}

bool xreal::operator > (const xreal y) {
    return (i > y.i || (i == y.i && f > y.f));
}

bool xreal::operator >= (const xreal y) {
    return (i > y.i || (i == y.i && f >= y.f));
}

real fmod2(xreal x, real y)	// limited function:  assumes y is a power
				//		      of 2 less than 1
{
    xfrac_t fx = x.get_f();
    xreal xy = (xreal)y;
    xfrac_t fy = xy.get_f();
    xfrac_t r = fx%fy;
    return r*TWO64I;
}

#else

real fmod2(xreal x, real y)
{
    return fmod(x, y);
}

#endif


//----------------------------------------------------------------------
//
// All the xreal I/O functions should appear here, for consistency.

// Input functions are derived from a single get_xreal function.

#ifdef USE_XREAL

#include <sstream>

local bool read_real_from_stream(istringstream& ss, real& x)
{
    ss.seekg(0, ios::beg);	// rewind
    return (ss >> x && (ss.eof() || ss.peek() <= ' '));
}

local bool read_int_from_stream(istringstream& ss, xint_t& i,
				real& x, bool& xset)
{
    // Read a (signed) decimal integer from the input stream.
    // Return true if something was successfully read.
    // Set xset = true if we read a real number and return
    // its value in x; otherwise we return the value in i.

    i = 0;
    xset = false;

    if (ss >> dec >> i) {

	// We have successfully read something...
	// Check to see if the next character is whitespace.
	// If not, then something unexpected occurred.

	if (!ss.eof() && ss.peek() > ' ' && ss.peek() != '+')

	    // Read stopped when it encountered an unexpected
	    // character.  Attempt to read a real number.

	    return (xset = read_real_from_stream(ss, x));

    } else {

	// Deal with error conditions.

	if (ss.bad()) return false;
	if (ss.fail()) ss.clear();

	// Read stopped when it encountered an unexpected character.

	if (ss.peek() != '+')

	    // Attempt to read a real number.

	    return (xset = read_real_from_stream(ss, x));
    }

    return true;
}

local bool read_uint_from_stream(istringstream& ss, xfrac_t& f)
{
    // Read an unsigned integer from the input stream.  No way
    // to tell the difference between dec 123 and hex 123 if the
    // leading 0x is not mandatory, so simply always assume that
    // the string is hex, with or without a 0x.  Return true iff
    // a valid value for f was successfully read.

    f = 0;

    // Skip leading whitespace and '+'...

    unsigned char c = 0;
    while (!ss.eof() && (c <= ' ' || c == '+')) ss >> c;
    if (ss.eof()) return true;
    ss.putback(c);

    if (ss >> hex >> f >> dec)

	// If we are not at whitespace or the end of the string,
	// then an error has occurred.

	return (ss.eof() || ss.peek() <= ' ');

    else

	return false;
}

local bool read_int_pair_from_stream(istringstream& ss,
				     xint_t& i, xfrac_t& f,
				     real& x, bool& xset)
{
    // Read an (int, unsigned int) pair from stream ss.
    // Excessively complicated code to cover the various input options.
    //
    //		int   uint (hex)
    //		int + uint (hex)
    //		int   0xuint
    //		int + 0xuint
    //		real

    static const char *func = "read_int_pair_from_stream";

    i = 0;
    f = 0;
    x = 0;
    xset = false;

    // First read a signed integer from the stream.

    if (read_int_from_stream(ss, i, x, xset)) {

	if (xset)

	    // Got a real number.

	    return  true;

	else {

	    // Got the first integer.  Now read a second (unsigned).

	    if (read_uint_from_stream(ss, f))
		return true;
	    else
		cerr << func << ": unable to read string (2) "
		     << ss.str() << endl;
	}

    } else

	cerr << func << ": unable to read string (1) " << ss.str() << endl;

    return false;
}

static bool print_xreal = true;

local inline xreal read_xreal(const char *str)
{
    // Extract an xreal from a string, with tests for various formats.
    // Should be able to understand
    //
    //		int   uint
    //		int + uint
    //		int   uint (hex)
    //		int + uint (hex)
    //		int   0xuint
    //		int + 0xuint
    //		real
    //
    // Not all options work with both the old and the new code...

#if 0

    // Old version uses strtoll, etc., which causes problems for configure.
    // This version doesn't understand "+" between the ints.

    char *sp, *ep;
    xint_t i = STRTOL(str, &sp, 10);		  // signed integer part
    xfrac_t f = STRTOUL(sp, &ep, 0); 		  // unsigned fractional part
						  // "0" here means that we
						  // can read hex or integer
    // PRC(f); PRL(ep);

    if (sp == ep) {				  // if we didn't get both
						  // of above,

	// Hmmm... most likely we have real input data.  Try just reading
	// a real number.  (Steve, 6/00)

	if (print_xreal) {
	    cerr << "read_xreal: error reading xreal input "
		 << "from string" << endl
		 << "    " << str << endl
		 << "Assuming real data." << endl << endl;
	    print_xreal = false;
	}

	return (xreal)strtod(str, NULL);
    }

    return xreal(i, f);

#else

    // New version is 100% C++, but (long) long!

    xint_t i;					  // signed integer part
    xfrac_t f;			 		  // unsigned fractional part
    real x;
    bool xset;

    istringstream ss(str);
    if (!read_int_pair_from_stream(ss, i, f, x, xset)) {
	if (print_xreal) {
	    cerr << "read_xreal: error reading xreal input "
		 << "from string" << endl
		 << "    " << str << endl << endl;
	    print_xreal = false;
	}
	return xreal(0, 0);
    }

    if (xset)
	return x;
    else
	return xreal(i, f);

#endif
}



// Note that the >> operator is not symmetric with <<...

istream & operator >> (istream & s, xreal & x)
{
#if 0

    // C++ input apparently doesn't understand "0xa123bcde" format for hex...

    // xint_t i;
    // xfrac_t f;
    // s >> i >> f;		// s >> x.i >> x.f fails; don't know why...
    // x = xreal(i, f);

    // Operator is not widely used; just use read_xreal().

    // Start by reading in two "words".  Easiest to use STL strings for this.

    string i, f;
    s >> i >> f;
    string str = i + " " + f;
    x = read_xreal(str.c_str());

#else

    // Hmmm.  Now we know how to make C++ do what we want...

    xint_t i;
    xfrac_t f;
    s >> i >> hex >> f >> dec;		// still needed?

    // Note that we currently don't check for input errors...
    // For a more thorough function to read an xreal from a
    // stream, see read_xreal().

    x = xreal(i, f);

#endif

    return s;
}

#endif

xreal get_xreal(char *str)
{
#if defined USE_XREAL

    // "True" xreal:

    return read_xreal(str);

#else

    // xreal is really just real:

    return (xreal)strtod(str, NULL);

#endif
}

xreal get_xreal_from_input_line(char * input_line)
{
    char *val = strchr(input_line, '=');
    if(val == NULL) return (xreal)0;
    val++;

#if defined USE_XREAL

    // "True" xreal:

    return read_xreal(val);

#else

    // xreal is really just real:

    return (xreal)strtod(val, NULL);

#endif
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// All xreal "print" functions are now just versions of a single xprint.

void xprint(xreal x,
	    ostream & s,	// default = cerr
	    bool newline)	// default = true
{
#ifdef USE_XREAL

#if 1

    // Simplest version:
    //
    // s << label << x.get_i() << " " << x.get_f() << endl;
    //
    // but now use hex (with leading 0x) for the fractional part.
    // Use 0x instead of showbase because g++ 2.95 doesn't understand.

    s << x.get_i() << " " << hex << "0x" << x.get_f() << dec;

#else

    // Older version: use C printf to do the job.

    xfrac_t f = x.get_f();
    char tmp[128];
    if (f == 0)
	sprintf(tmp, "0");		// handy
    else
	sprintf(tmp, "%#16.16llx", f);

    s << x.get_i() << " " << tmp;

#endif

#else

    s << x;

#endif

    if (newline) s << endl;
}

#ifdef USE_XREAL

// Member function xprint is rarely used, but convenient to retain it.
// Simply define it in terms of the non-member xprint function.

void xreal::print(ostream& s,
		  bool newline)	// default = false
{
    xprint(*this, s, newline);
}

// Xreal version of put_real_number is just xprint with a label.

void put_real_number(ostream & s, char * label, xreal x)
{
    s << label;
    xprint(x, s);
}
#endif

void identify_xreal(ostream &s)	// default = cerr
{
    s << "xreal data type is ";
#ifdef USE_XREAL
    s << "extended-precision real" << endl;
#else
    s << "real" << endl;
#endif
}

//----------------------------------------------------------------------


#else

main()
{
    identify_xreal();

#ifdef USE_XREAL

    cerr.precision(HIGH_PRECISION);

    while (1) {

	char s[256];
	cout << endl << "enter x.i x.f or real: " << flush;

	if (!cin.getline(s, sizeof(s)-1)) break;

	xreal x = get_xreal(s);
	put_real_number(cerr, "x = ", x);

	cerr << "real "; PRC(x); cerr << "xreal: "; x.print(); cerr << endl;

	xreal y, z;
	z = (y = x + 1.0);
	cerr << "real "; PRC(y); cerr << "xreal: "; y.print(); cerr << endl;
	cerr << "real "; PRC(z); cerr << "xreal: "; z.print(); cerr << endl;

	real dx = 1.e-18;
	y = x + dx;
	cerr << "real "; PRC(y); cerr << "xreal: "; y.print(); cerr << endl;

	real dy = y - x;
	PRL(dy);
    }

#endif
}

#endif
