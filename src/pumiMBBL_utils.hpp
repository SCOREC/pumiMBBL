#ifndef pumiMBBL_utils_hpp
#define pumiMBBL_utils_hpp

namespace pumi{

/**
 * @brief Copies a derived class to device and returns pointer to it.
 *
 * [COPIED FROM HPIC2 UTILS -- NECESSARY FOR SUBMESH-LEVEL POLYMORPHIC APIs ]
 * Necessary for using virtual functions on GPUs.
 * This mallocs data for the Derived class on Device,
 * then returns a Base pointer to that data.
 * This allows us to call virtual functions from the device.
 *
 * @note Since this is dynamically malloc'd,
 *       the data must be free'd before conclusion.
 *
 * @param[in] to_copy Derived class to be copied to device.
 *
 * @return Base pointer to copied object on device.
 */
template<typename Base, typename Derived, typename Device=Kokkos::DefaultExecutionSpace>
Base* copyForDevice(const Derived & to_copy) {
    auto * p = static_cast<Base *>(Kokkos::kokkos_malloc<typename Device::memory_space>(sizeof(Derived)));
    Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0, 1), KOKKOS_LAMBDA (const int) {
        new (p) Derived(to_copy);
    });
    Kokkos::fence();
    return p;
}

/**
 * @brief Frees T object from device.
 *
 * After having dynamically copied an object to device,
 * it must be manually free'd.
 *
 * @param[in] ptr Pointer to object to be deleted.
 */
template<typename T, typename Device=Kokkos::DefaultExecutionSpace>
void deleteOnDevice(T *ptr) {
    // It may not be necessary to to kokkos_free from within a parallel_for.
    Kokkos::kokkos_free<typename Device::memory_space>(ptr);
    // Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0, 1),
    //                     KOKKOS_LAMBDA (const int) {
    //     Kokkos::kokkos_free<typename Device::memory_space>(ptr);
    // });
}

/**
 * @brief Convenience class for containing pointers in views.
 */
template<typename T>
class DevicePointer {
    public:
        /// Default constructor
        KOKKOS_INLINE_FUNCTION
        DevicePointer() : ptr_(nullptr) {}

        /// Default destructor
        ~DevicePointer() = default;

        /// Constructor saves existing pointer
        KOKKOS_INLINE_FUNCTION
        DevicePointer(T *ptr) : ptr_(ptr) {}

        /// Return pointer
        KOKKOS_INLINE_FUNCTION
        T * operator() () { return ptr_; }

    private:
        T *ptr_;
};



}
#endif
