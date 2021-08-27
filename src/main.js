const core = require('@actions/core')
const doInit = require('./init.js')
const doEval = require('./eval.js')

try {
  const task = core.getInput('task')
  if (task === 'download') {
    doInit().catch(error => core.setFailed(error.message))
  } else if (task === 'eval') {
    doEval().catch(error => core.setFailed(error.message))
  } else {
    core.setFailed(`Invalid task ${task}`)
  }
} catch (error) {
  core.setFailed(error.message)
}
