def cisynergy(distribution, return_union_information = False):
    """This function receives a variable of type dit.Distribution and computes the union information of the sources I_\cap (Y_1, ..., Y_n -> T)
    for n=2 or n=3. On the input distribution, the target T is assumed to be the last variable. It requires the dit package to run. If return_synergy is True,
    this function also returns synergy S, as defined by S = I(Y;T) - I_\cap (Y_1, ..., Y_n -> T).
    
    One must first check if any of the sources Y_i is a deterministic function of other source Y_j, and if so, source Y_i must be removed manually from
    the dit.Distribution before computation. This can be easily checked since Y_i is a deterministic function of Y_j iff H(Y_i|Y_j) = 0.
    """
    import dit
    
    d = distribution.copy()
    
    if d.outcome_length() == 3:
        n_sources = 2
    elif d.outcome_length() == 4:
        n_sources = 3
    else:
        return None
    
    if n_sources == 2:

        d.set_rv_names('XYT')

        Tsupport = list(d.marginal('T').outcomes)
        Xsupport = list(d.marginal('X').outcomes)
        Ysupport = list(d.marginal('Y').outcomes)

        pT = d.marginal('T')
        pXT = d.marginal('XT')
        pYT = d.marginal('YT')

        Qsupport = []

        for i in Tsupport:
            for j in Xsupport:
                for k in Ysupport:
                    Qsupport.append(j+k+i)

        Qprobabilities = []

        for i in Qsupport:
            X = i[0]
            Y = i[1]
            T = i[2]
            try:
                prob = pT[T] * pXT[X+T]/pT[T] * pYT[Y+T]/pT[T]
            except:
                prob = 0
            Qprobabilities.append(prob)

        distribuicaoQ = dit.Distribution(Qsupport, Qprobabilities)

        distribuicaoQ.set_rv_names('XYT')

        union_information = min(dit.shannon.mutual_information(distribuicaoQ, ['T'], ['X', 'Y']), dit.shannon.mutual_information(d, ['T'], ['X', 'Y']))
        
        synergy = dit.shannon.mutual_information(d, ['T'], ['X', 'Y']) - union_information
        
        if return_union_information:
            return [synergy, union_information]
        else:
            return synergy
        
    elif n_sources == 3:
        
        d.set_rv_names('XYZT')

        Tsupport = list(d.marginal('T').outcomes)
        Xsupport = list(d.marginal('X').outcomes)
        Ysupport = list(d.marginal('Y').outcomes)
        suporteZ = list(d.marginal('Z').outcomes)

        pT = d.marginal('T')
        pXT = d.marginal('XT')
        pYT = d.marginal('YT')
        pZT = d.marginal('ZT')

        Qsupport = []

        for i in Tsupport:
            for j in Xsupport:
                for k in Ysupport:
                    for l in suporteZ:
                        Qsupport.append(j+k+l+i)

        Qprobabilities = []

        for i in Qsupport:
            X = i[0]
            Y = i[1]
            Z = i[2]
            T = i[3]
            try:
                prob = pT[T] * pXT[X+T]/pT[T] * pYT[Y+T]/pT[T] * pZT[Z+T]/pT[T]
            except:
                prob = 0
            Qprobabilities.append(prob)

        distribuicaoQ = dit.Distribution(Qsupport, Qprobabilities)

        distribuicaoQ.set_rv_names('XYZT')

        union_information = min(dit.shannon.mutual_information(distribuicaoQ, ['T'], ['X', 'Y', 'Z']), dit.shannon.mutual_information(d, ['T'], ['X', 'Y', 'Z']))

        synergy = dit.shannon.mutual_information(d, ['T'], ['X', 'Y', 'Z']) - union_information

        if return_union_information:
            return [synergy, union_information]
        else:
            return synergy