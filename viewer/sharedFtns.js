
function unique (L) {
  // get unique values in list
  var toRet = new Set(L)
  return Array.from(toRet)
}

if (typeof module !== 'undefined') {
  module.exports = { unique }
}
