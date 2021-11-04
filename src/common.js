const core = require('@actions/core')
const exec = require('@actions/exec')
const yaml = require('js-yaml')
const fs = require('fs')
const path = require('path')

const benchmarkName = core.getInput('benchmark_name')


const benchmarkOutdir = `benchmark-data/${benchmarkName}`
const conda = `${benchmarkOutdir}/mamba/bin/conda`
const mamba = `${benchmarkOutdir}/mamba/bin/mamba`

async function _exec (cmd) {
  await exec.exec('bash', ['-l', '-c', `${cmd}`])
}

let meta = path.join(__dirname, `../benchmarks/${benchmarkName}/meta.yaml`)
meta = yaml.load(fs.readFileSync(meta, 'utf-8'))

const common = {
  exec: async function (cmd) {
    await _exec(cmd)
  },
  initMamba: async function () {
    await _exec(`curl -L --insecure https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh > mambaforge.sh; bash mambaforge.sh -b -p ${benchmarkOutdir}/mamba`)
    await _exec(`${benchmarkOutdir}/mamba/bin/conda init bash`)
  },
  getEnvActivate: function (envpath, name) {
    if (fs.existsSync(envpath)) {
      return `${conda} activate ${name}`
    }
    return ''
  },
  initEnv: async function (envpath, name) {
    if (fs.existsSync(envpath)) {
      await _exec(`source ${benchmarkOutdir}/mamba/etc/profile.d; echo $PATH; ${conda} activate base; ${mamba} env create --yes -n ${name} -f ${envpath}`)
    }
  },
  getBenchmarkName: function () {
    return benchmarkName
  },
  getBenchmarkFile: function (filename) {
    return path.join(__dirname, `../benchmarks/${benchmarkName}/${filename}`)
  },
  getBenchmarkOutdir: function () {
    return benchmarkOutdir
  },
  getOutpathEnvvars: function () {
    return Object.entries(meta.paths).map(function ([key, path]) { return `${key}="${benchmarkOutdir}/${path}"` }).join(' ')
  },
  getOutpath: function (name) {
    return `${benchmarkOutdir}/${meta.paths[name]}`
  }
}

module.exports = common
