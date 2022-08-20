const fs = require('fs')
const path = require('path')

const iter = parseInt(process.argv[2])
const tqpPath = process.argv[3]
const allqpPath = process.argv[4]
const outPath = process.argv[5]

const tqp = JSON.parse(fs.readFileSync(path.resolve(tqpPath), 'utf-8'))
const allqp = JSON.parse(fs.readFileSync(path.resolve(allqpPath), 'utf-8'))

const qpFtns = require('./quiverFtns.js')

var result = qpFtns.potentialRandomSearch(tqp, allqp.length, allqp, 4, 1, 30, iter)

fs.writeFileSync(path.resolve(outPath), JSON.stringify(result))
