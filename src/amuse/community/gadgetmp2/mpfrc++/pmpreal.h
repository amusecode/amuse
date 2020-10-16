/*    
    Copyright (c) 2020 Thomas Schano

    Contributors:

    Licensing:
    (A) PMPREAL is under GNU General Public License ("GPL").

    (B) Non-free licenses may also be purchased from the author, for users who
        do not want their programs protected by the GPL.

        The non-free licenses are for users that wish to use MPFR C++ in
        their products but are unwilling to release their software
        under the GPL (which would require them to release source code
        and allow free redistribution).

        Such users can purchase an unlimited-use license from the author.
        Contact us for more details.

    GNU General Public License ("GPL") copyright permissions statement:
    **************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PMPREAL_H_INCLUDED
#define PMPREAL_H_INCLUDED
#include "mpreal.h"

/*
 * This file needs a modified mpreal.h!
 * The class mpreal needs a public friend reference to pmpreal
 * ->add following line to public section:  friend class pmpreal;
 *
 */

namespace mpfr {
class pmpreal : public mpreal
{
private:
    inline void* operator new[](size_t size);
    void* operator new(size_t size)= delete;
    void  operator delete(void*)= delete;
    void  set_prec() = delete;
    mpreal& setPrecision() = delete;
public:
    inline pmpreal& operator=(const mpreal& v);
    static inline size_t get_needed_mem(size_t size);
    inline size_t get_transfer_buff_size();
    static inline size_t get_needed_mem_single(mpfr_prec_t prec);
    static inline pmpreal* generate_buffer(size_t size);
    static inline pmpreal* place_buffer(size_t size, void* buff_start, mpfr_prec_t prec);
    static inline void place_pmpreal(void* pmpreal, mpfr_prec_t prec);
    inline size_t size();
    static inline char* get_transfer_buff(pmpreal* p);
    inline char* get_transfer_buff();
    static inline size_t get_buff_size(pmpreal* p);
    static inline void re_org_buff(pmpreal* p);
    static inline pmpreal* re_org_transfer_buff(void* p);
    static inline void re_org_pmpreal(pmpreal* p);
    inline void* next_pointer();
    inline ~pmpreal();
};

size_t pmpreal::get_transfer_buff_size()
{
    return get_needed_mem(size());
}

pmpreal* pmpreal::re_org_transfer_buff(void* p)
{
    pmpreal* retval=(pmpreal*)((size_t)p+sizeof(size_t));
    re_org_buff(retval);
    return retval;
}
void pmpreal::re_org_buff(pmpreal* p)
{
    size_t size=p->size();
    size_t headers = size*sizeof(pmpreal);
    size_t mantissa = mpfr_custom_get_size(p[0].getPrecision());
    size_t new_p = (size_t)p + headers;
    void* new_p_p;
    for (size_t i=0;i<size;i++)
    {
        new_p_p=(void*)new_p;
        (mpfr_custom_move)(p[i].mp,new_p_p);
        new_p+=mantissa;
    };

}

void pmpreal::re_org_pmpreal(pmpreal* p)
{
    size_t header = sizeof(pmpreal);
    size_t mantissa = mpfr_custom_get_size(p->getPrecision());
    size_t new_p = (size_t)p + header;
    void* new_p_p;
        new_p_p=(void*)new_p;
        (mpfr_custom_move)(p->mp,new_p_p);
}

size_t pmpreal::get_buff_size(pmpreal* p)
{
    size_t mantissa=mpfr_custom_get_size(mpfr_get_prec(p[0].mp));
    return (sizeof(pmpreal)+mantissa)*p->size() + sizeof(size_t);
}

void pmpreal::place_pmpreal(void* pointer, mpfr_prec_t prec=0)
{
    mpfr_prec_t preccheck;
    if (prec == 0)
        preccheck=mpreal::get_default_prec();
    else
        preccheck=prec;
    size_t header = sizeof(pmpreal);
    void* new_p_p = (void*)((size_t)pointer +   header);
    pmpreal* tmp = (pmpreal*) pointer;
     (mpfr_custom_init_set)(tmp->mp, MPFR_ZERO_KIND, 1,preccheck,new_p_p);
    return;
}



pmpreal* pmpreal::place_buffer(size_t size, void* buff_start, mpfr_prec_t prec=0)
{
    mpfr_prec_t preccheck;
    if (prec == 0)
        preccheck=mpreal::get_default_prec();
    else
        preccheck=prec;
    size_t* arr_size=  (size_t*) buff_start;
    *arr_size = size;
    pmpreal* retval = (pmpreal*)(arr_size+1);
    size_t headers = size*sizeof(pmpreal);
    size_t mantissa = mpfr_custom_get_size(preccheck);
    size_t pointe = (size_t)retval;
    size_t new_p = pointe + headers;
    void* new_p_p;
    for (size_t i=0;i<size;i++)
    {
        new_p_p=(void*)new_p;
        (mpfr_custom_init_set)(retval[i].mp, MPFR_ZERO_KIND, 1,preccheck,new_p_p);
        new_p+=mantissa;
       };
    return retval;
}

pmpreal* pmpreal::generate_buffer(size_t size)
{
    mpfr_prec_t preccheck=mpreal::get_default_prec();
    pmpreal* retval = new pmpreal[size];
    if(retval[0].getPrecision()==preccheck) // other thread may have changed the default prec
    {
       size_t headers = size*sizeof(pmpreal);
       size_t mantissa = mpfr_custom_get_size(preccheck);
       size_t pointe = (size_t)retval;
       size_t new_p = pointe + headers;
       void* new_p_p;
       for (size_t i=0;i<size;i++)
       {
            new_p_p=(void*)new_p;
            mpfr_clear(retval[i].mp);
            (mpfr_custom_init_set)(retval[i].mp, MPFR_ZERO_KIND, 1,preccheck,new_p_p);
            new_p+=mantissa;
       };
    }
    else{retval=NULL;}
    return retval;
}

pmpreal::~pmpreal()
{
    mpfr_init(mp);
}

size_t pmpreal::size()
{
    size_t *p = ((size_t*)this)-1;
    return (size_t)(*p);
}


void* pmpreal::next_pointer()
{
    pmpreal *tmp= (pmpreal*)this;
    size_t retval=(sizeof(pmpreal)+mpfr_custom_get_size(tmp[0].getPrecision()))*tmp->size();
    return (void*)((size_t)tmp+retval);
}

char* pmpreal::get_transfer_buff()
{
    pmpreal* p= (pmpreal*)this;
    char* retval =(char*)(((size_t*)p)-1);
    return retval;
}

char* pmpreal::get_transfer_buff(pmpreal* p)
{
    char* retval =(char*)(((size_t*)p)-1);
    return retval;
}

pmpreal& pmpreal::operator=(const mpreal& v)
{
    if (this != &v)
    {
        mpfr_set(mpreal::mp, v.mp, mpreal::get_default_rnd());
    }
    return *this;
}

void * pmpreal::operator new[](size_t size)
{
    size_t counts=(size-sizeof(size_t))/sizeof(pmpreal) ;
	size_t limbs=mpfr_custom_get_size(mpreal::get_default_prec())*counts;
	void * p = malloc(size+limbs);
	return p;
}


size_t pmpreal::get_needed_mem(size_t size)
{
    size_t limbs=mpfr_custom_get_size(mpreal::get_default_prec())*size;
    size_t header= size* sizeof(pmpreal) + sizeof(size_t);
    return header + limbs;
}

size_t pmpreal::get_needed_mem_single(mpfr_prec_t prec=0)
{
    mpfr_prec_t preccheck;
    if (prec == 0)
        preccheck=mpreal::get_default_prec();
    else
        preccheck=prec;

    size_t limbs=mpfr_custom_get_size(preccheck);
    size_t header= sizeof(pmpreal);
    return header + limbs;
}
}
#endif // PMPREAL_H_INCLUDED
