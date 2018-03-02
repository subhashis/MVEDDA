
#ifndef SHARED_ARRAY_H_
#define SHARED_ARRAY_H_

#include <memory>

namespace edda{

///
/// smart pointer for arrays
///
template <class T>
class shared_ary : public std::shared_ptr<T>
{
    size_t n;
    class DeleteArray
    {
    public:
        void operator () (T* d) const
        {
          if (d)
            delete [] d;
        }
    };

public:
    shared_ary(): std::shared_ptr<T>(NULL, DeleteArray()), n(0) {}
    shared_ary(T *p, size_t n_): std::shared_ptr<T>(p, DeleteArray()), n(n_) { }
    shared_ary(size_t n_): shared_ary(new T[n_], n_) {}

    T & operator[] (int i) const { //std::cout << n << ',' << i;
                                   assert(i>=0 && i<n); return this->get()[i]; }
    size_t getLength() const {return n;}
    void swap(shared_ary<T> &ary) {
      std::shared_ptr<T>::swap(ary); std::swap(n, ary.n);
    }
};

} // namespace edda

#endif // SHARED_ARRAY_H_
