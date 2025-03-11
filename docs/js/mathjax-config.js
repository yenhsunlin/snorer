window.MathJax = {
  tex: {
    packages: ['base', 'ams'], // 確保 AMS 環境可用
    tags: 'ams',  // 讓標號自動對齊
    macros: {
      subeq: ['{\\begin{subequations} #1 \\end{subequations}}', 1]
    }
  },
  loader: {load: ['[tex]/ams', '[tex]/autoload']},
};