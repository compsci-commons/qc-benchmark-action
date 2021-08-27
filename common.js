const exec = require('@actions/exec');



var benchmarkName = core.getInput('benchmark-name');

var benchmarkOutdir = `benchmark-data/${benchmarkName}`;
var bashrc = `${benchmarkOutdir}/.bashrc`;
var micromamba = `${benchmarkOutdir}/bin/micromamba`;
var bashPrefix = `source ${bashrc}; `;

fs.closeSync(fs.openSync(bashrc, 'w'))

function _exec(cmd) {
    return exec('bash', ['-l', '-c', `source ${bashrc}; ${command}`]);
}

var common = {
    exec: function(cmd) {
        return _exec(cmd);
    },
    initMicromamba: function() {
        _exec(`curl -L https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj ${micromamba}`);
        _exec(`${micromamba} shell init -s bash -p ${benchmarkOutdir}/tmp/micromamba --rc-file ${bashrc}`);
    },
    initEnv: function(envpath, name) {
        let prefix = '';
        if(fs.existsSync(envpath)) {
            _exec(`${micromamba} create -n ${name} -f ${envpath}`);
            return `micromamba activate ${name}`;
        }
        return '';
    },
    getBenchmarkName: function() {
        return benchmarkName;
    },
    getBenchmarkFile: function(filename) {
        return `benchmarks/${benchmarkName}/${filename}`;
    },
    getBenchmarkOutdir: function() {
        return benchmarkOutdir;
    }
}

module.exports = common;