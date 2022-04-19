function displayMinimalFor(results, exchangeNum, idx=0) {
    potential.clear()
    potential.add(results.minimalResultsByExchangeNum[exchangeNum][idx].map(function(term) {
        return {id: term[1], coef: term[0].toString()}
    }))
}
/*
function shrinkPotential(qp, exchangeNum, potential) {
    var qpt = deepCopy(qp)
    qpt.potential = potential

    var canShrink = true;
    while (canShrink) {
        canShrink = false

        for (var i = 0; i < potential.length; i++) {
            var potentialWithoutI = potential.filter((t, idx) => idx != i)

            var qpt2 = deepCopy(qpt)
            qpt2.potential = potentialWithoutI

            var mut = getAllMutationsForQP(qpt2);

            if (mut.quivers.length === exchangeNum) {
                qpt = qpt2;
                canShrink = true;
                break;
            }
        }
    }
    return qpt.potential
}
*/