#ifndef _vec_hpp_INCLUDED
#define _vec_hpp_INCLUDED

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstddef>

template<class T>
class vect {
    T* data;
    static inline int  imax   (int x, int y) { int mask = (y-x) >> (sizeof(int)*8-1); return (x&mask) + (y&(~mask)); }
    static inline void nextCap(int &cap) { cap += ((cap >> 1) + 2) & ~1; }

    public:
        int sz, cap;
        vect()                   :  data(NULL), sz(0), cap(0) {}
        explicit vect(int size)  :  data(NULL), sz(0), cap(0) { growTo(size); }
        ~vect()                                               { clear(true); }

        operator T*      (void)         { return data; }
        
        int     size     (void) const   { return sz;   }
        int     capacity (void) const   { return cap;  }

        void    setsize  (int v)        { sz = v;} 
        void    push   (const T& elem)  { if (sz == cap) capInc(sz + 1); data[sz++] = elem; }
        void    push_  (const T& elem)  { assert(sz < cap); data[sz++] = elem; }
        void    pop    (void)           { assert(sz > 0), sz--, data[sz].~T(); } 
        void    copyTo (vect<T>& copy)   { copy.clear(); copy.growTo(sz); for (int i = 0; i < sz; i++) copy[i] = data[i]; }
        
        void    growTo   (int size);
        void    clear    (bool dealloc = false);
        void    capInc   (int to_cap);


        T&       operator [] (int index)       { return data[index]; }
        const T& operator [] (int index) const { return data[index]; }

        T&       last        (void)            { return data[sz - 1]; }
        const T& last        (void)      const { return data[sz - 1]; }
                    
};


class OutOfMemoryException{};

template<class T>
void vect<T>::clear(bool dealloc) {
    if (data != NULL) {
        sz = 0;
        if (dealloc) free(data), data = NULL, cap = 0;
    }
}

template<class T>
void vect<T>::capInc(int to_cap) {
    if (cap >= to_cap) return;
    int add = imax((to_cap - cap + 1) & ~1, ((cap >> 1) + 2) & ~1); 
    if (add > __INT_MAX__ - cap || ((data = (T*)::realloc(data, (cap += add) * sizeof(T))) == NULL) && errno == ENOMEM)
        throw OutOfMemoryException();
}

template<class T>
void vect<T>::growTo(int size) {
    if (sz >= size) return;
    capInc(size);
    for (int i = 0; i < sz; i++) new (&data[i]) T();
    sz = size;
}

#endif