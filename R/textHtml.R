tableau_html <- function(data, data2) {
  effets <- c("\\(\\zeta\\)", "\\(\\delta\\)")
  
  headers <- c(
    "Effets", "Méthode", "Corrélation estimée",
    "Valeur estimée", "Écart-type<sup>1</sup>", "95% IC<sup>1</sup>"
  )
  
  headers_html <- sapply(headers, function(h) {
    words <- unlist(strsplit(h, " "))
    if (length(words) > 1) {
      paste(words[1], "<br>", paste(words[-1], collapse = " "), sep = "")
    } else {
      h
    }
  })
  
  html <- '
  <style>
    table.custom-table {
      margin: auto;
      border-collapse: collapse;
      font-family: Arial, sans-serif;
      font-size: 14px;
      width: 80%;
      border: none;
    }
    .custom-table th,
    .custom-table td {
      text-align: center;
      vertical-align: middle;
      padding: 8px 12px;
      border: none;
    }
    .custom-table thead th {
      border-bottom: 2px solid #999999;
      background-color: #f2f2f2;
    }
    .custom-table tbody tr:nth-child(even) td {
      background-color: #f9f9f9;
    }
    /* Enlever les lignes horizontales sauf celle avant les résultats de data2 */
    .custom-table tbody tr:not(:first-child):not(:nth-child(3)) td {
      border-bottom: none;
    }
    .custom-table tbody tr:nth-child(3) td {
      border-bottom: 1px solid #dddddd; /* Ligne horizontale avant les résultats de data2 */
    }
    .custom-table tbody tr:first-child td {
      border-bottom: none;
    }
    .custom-table tbody tr:last-child td {
      border-bottom: none;
    }
  </style>

  <table class="custom-table">
    <thead>
      <tr>'
  
  for (h in headers_html) {
    html <- paste0(html, '<th>', h, '</th>')
  }
  
  html <- paste0(html, '</tr>
    </thead>
    <tbody>')
  
  # Première méthode : CC
  for (i in seq_along(effets)) {
    html <- paste0(html, '<tr><td>', effets[i], '</td>')
    if (i == 1) {
      html <- paste0(html,
                     '<td rowspan="2" style="vertical-align:middle;">\\(\\texttt{CC}\\)</td>')
    }
    html <- paste0(html, '<td>', "", '</td>')
    html <- paste0(html, '<td>', data[i, "Valeur"], '</td>')
    html <- paste0(html, '<td>', data[i, "ET"], '</td>')
    html <- paste0(html, '<td>', data[i, "IC"], '</td></tr>')
  }
  
  # Méthodes CNCm et CNCr
  effets2 <- c("\\(\\zeta\\)", "\\(\\delta\\)")
  
  for (effet in effets2) {
    sous_data <- data2[data2$Effet == effet, ]
    n <- nrow(sous_data)
    
    for (i in seq_len(n)) {
      ligne <- sous_data[i, ]
      
      if (i == 1) {  # Première ligne => CNCm
        html <- paste0(html, "<tr><td>", effet, "</td>")
        html <- paste0(html, '<td >\\(\\texttt{CNCm}\\)</td>')
      } else {  # Lignes suivantes => CNCr
        html <- paste0(html, "<tr><td></td>")
        html <- paste0(html, '<td rowspan="1" style="vertical-align:middle;">\\(\\texttt{CNCr}\\)</td>')
      }
      html <- paste0(html,
                     "<td>", ligne$rho, "</td>",
                     "<td>", ligne$Valeur, "</td>",
                     "<td>", ligne$ET, "</td>",
                     "<td>", ligne$IC, "</td></tr>")
    }
  }
  
  html <- paste0(html, '</tbody></table>')
  # Ajouter une ligne fine avant "Bootstrap"
  html <- paste0(html, '<hr style="border: 0; border-top: 1px solid #ddd;">')
  
  # Le texte "Bootstrap" à gauche et une ligne après
  html <- paste0(html, '<p style="text-align:left; margin-left: 20px;"><sup>1</sup> Bootstrap</p>')
  html <- paste0(html, '<hr style="border: 0; border-top: 1px solid #ddd;">')
  
  # Ajout de MathJax pour rendre le LaTeX
  html <- paste0(html, '
  <script type="text/javascript" async
    src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
  </script>
  <script type="text/javascript">
    MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
  </script>')
  
  return(HTML(html))
}