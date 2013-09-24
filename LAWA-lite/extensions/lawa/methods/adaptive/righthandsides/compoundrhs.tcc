namespace lawa {

template <typename T, typename Index, typename FirstRhs, typename SecondRhs, typename ThirdRhs,
          typename FourthRhs, typename FifthRhs>
CompoundRhs<T,Index,FirstRhs,SecondRhs,ThirdRhs,FourthRhs,FifthRhs>::
CompoundRhs(FirstRhs  &_firstRhs, SecondRhs &_secondRhs)
:   numOfRhs(2), firstRhs(_firstRhs), secondRhs(_secondRhs), thirdRhs(_secondRhs),
    fourthRhs(_secondRhs), fifthRhs(_secondRhs)
{

}

template <typename T, typename Index, typename FirstRhs, typename SecondRhs, typename ThirdRhs,
          typename FourthRhs, typename FifthRhs>
CompoundRhs<T,Index,FirstRhs,SecondRhs,ThirdRhs,FourthRhs,FifthRhs>::
CompoundRhs(FirstRhs  &_firstRhs, SecondRhs &_secondRhs, ThirdRhs &_thirdRhs)
:   numOfRhs(3), firstRhs(_firstRhs), secondRhs(_secondRhs), thirdRhs(_thirdRhs),
    fourthRhs(_secondRhs), fifthRhs(_secondRhs)
{

}

template <typename T, typename Index, typename FirstRhs, typename SecondRhs, typename ThirdRhs,
          typename FourthRhs, typename FifthRhs>
T
CompoundRhs<T,Index,FirstRhs,SecondRhs,ThirdRhs,FourthRhs,FifthRhs>::operator()(const Index &index)
{
    T val = 0.;
    switch (numOfRhs)
    {
        case 2:
            val = firstRhs(index);
            val+= secondRhs(index);
            return val;
            break;
        case 3:
            val = firstRhs(index);
            val+= secondRhs(index);
            val+= thirdRhs(index);
            return val;
            break;
        default:
            std::cerr << "CompoundRhs not yet implemented for " << numOfRhs
                      << " operators. Exit." << std::endl;
            exit(1);
    }
}


}   // namespace lawa

