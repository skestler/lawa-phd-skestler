namespace flens {

template <typename T>
template <typename IndexType>
T *
DefaultAllocator<T>::allocate(IndexType length, IndexType firstIndex)
{
    return static_cast<T *>(malloc(length*sizeof(T))) - firstIndex;
}

template <typename T>
template <typename IndexType>
void
DefaultAllocator<T>::release(T *&data, IndexType firstIndex)
{
    free(data+firstIndex);
    data=0;
}

} // namespace flens
