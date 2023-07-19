const n = 4;
let matrix = [];
let savedMatrix = [];
let savedAlgebra = [];

function initMatrix(n) {
    matrix = new Array(n).fill().map(() => new Array(n).fill('0'));
}

function updateMatrix() {
    const matrixDiv = document.getElementById('matrix');
    matrixDiv.innerHTML = '';
    matrixDiv.style.gridTemplateColumns = `repeat(${n}, 40px)`;
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            const input = document.createElement('input');
            input.value = matrix[i][j];
            input.onchange = function() {
                matrix[i][j] = this.value;
            }
            matrixDiv.appendChild(input);
        }
    }
}

document.getElementById('increase').onclick = function() {
    n++;
    initMatrix(n);
    updateMatrix();
}

document.getElementById('decrease').onclick = function() {
    if (n > 1) {
        n--;
        initMatrix(n);
        updateMatrix();
    }
}

document.getElementById('addMatrix').onclick = function() {
    savedMatrix.push(matrix);
    const matrixPreview = document.getElementById('matrixPreview');
    
    const formattedMatrix = '[' + matrix.map(row => row.join(', ')).join(';') + ']';
    matrixPreview.value += formattedMatrix + '\n';

    initMatrix(n);
    updateMatrix();
}

document.getElementById('addAlgebra').onclick = function() {
    const matrixPreview = document.getElementById('matrixPreview');
    matrixPreview.value += '@\n';
    initMatrix(n);
    updateMatrix();
}

document.getElementById('clearMatrix').onclick = function() {
    initMatrix(n);
    updateMatrix();
}

document.getElementById('saveContinue').onclick = function() {
    alert("Saved to MatrixInput.txt");
}
document.getElementById('mode').addEventListener('change', function(e) {
    document.body.classList.toggle('night-mode', e.target.checked);
});

initMatrix(n);
updateMatrix();


