function analyzePotentials(results) {
    var ret = []
    for(var num in results.resultsByExchangeNum) {
        var sorted = results.resultsByExchangeNum[num].sort((a, b) => a.length - b.length)
        var filtered = sorted.filter(i => i.length === sorted[0].length)

        ret.push([num, filtered])
    }
    return ret;
}

function displayMinimalFor(results, exchangeNum, idx=0) {
    var sorted = results.resultsByExchangeNum[exchangeNum].sort((a, b) => a.length - b.length)
    var filtered = sorted.filter(i => i.length === sorted[0].length)

    potential.clear()
    potential.add(filtered[idx].map(function(term) {
        return {id: term[1], coef: term[0].toString()}
    }))
}