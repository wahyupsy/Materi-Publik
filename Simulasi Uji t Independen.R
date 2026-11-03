library(shiny)
library(ggplot2)

hitung_uji <- function(g1, g2, var_equal, alpha) {
  tt <- t.test(g1, g2, var.equal = var_equal, conf.level = 1 - alpha)
  n1 <- length(g1)
  n2 <- length(g2)
  selisih <- mean(g1) - mean(g2)
  sp <- sqrt(((n1 - 1) * var(g1) + (n2 - 1) * var(g2)) / (n1 + n2 - 2))
  d <- selisih / sp
  keputusan <- if (tt$p.value < alpha) "H0 ditolak" else "H0 tidak ditolak"

  tabel <- data.frame(
    Statistik = c("Selisih mean (K1 - K2)", "t", "df", "p-value",
                  paste0("Batas bawah CI ", round((1 - alpha) * 100), "% (selisih)"),
                  "Batas atas CI (selisih)", "Cohen's d", "Keputusan"),
    Nilai = c(sprintf("%.4f", selisih),
              sprintf("%.4f", unname(tt$statistic)),
              sprintf("%.4f", unname(tt$parameter)),
              sprintf("%.4f", tt$p.value),
              sprintf("%.4f", tt$conf.int[1]),
              sprintf("%.4f", tt$conf.int[2]),
              sprintf("%.4f", d),
              keputusan)
  )

  ringkas <- data.frame(
    n1 = n1, n2 = n2,
    `Selisih mean` = selisih,
    t = unname(tt$statistic),
    df = unname(tt$parameter),
    `p-value` = tt$p.value,
    `Cohen's d` = d,
    Keputusan = keputusan,
    check.names = FALSE
  )
  list(tabel = tabel, ringkas = ringkas)
}

desk <- function(x, nama, alpha, mu, sigma) {
  n <- length(x)
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  tcrit <- qt(1 - alpha / 2, df = n - 1)
  data.frame(
    Kelompok        = nama,
    `Mean populasi` = mu,
    `Mean sampel`   = m,
    `SD populasi`   = sigma,
    `SD sampel`     = sd(x),
    n               = n,
    SE              = se,
    `CI bawah`      = m - tcrit * se,
    `CI atas`       = m + tcrit * se,
    Median          = median(x),
    Min             = min(x),
    Maks            = max(x),
    check.names     = FALSE
  )
}

uji_normalitas <- function(g1, g2, nama_riset, alpha) {
  s1 <- shapiro.test(g1)
  s2 <- shapiro.test(g2)
  p <- c(s1$p.value, s2$p.value)
  data.frame(
    Riset = nama_riset,
    Kelompok = c("Kelompok 1", "Kelompok 2"),
    W = c(unname(s1$statistic), unname(s2$statistic)),
    `p-value` = p,
    Keputusan = ifelse(p < alpha, "Tidak normal", "Normal"),
    check.names = FALSE
  )
}

uji_levene <- function(g1, g2, nama_riset, alpha) {
  x <- c(g1, g2)
  grp <- factor(rep(c("K1", "K2"), c(length(g1), length(g2))))
  ctr <- tapply(x, grp, median)
  z <- abs(x - ctr[as.character(grp)])
  m <- tapply(z, grp, mean)
  grand <- mean(z)
  ssb <- sum(tapply(z, grp, length) * (m - grand)^2)
  ssw <- sum((z - m[as.character(grp)])^2)
  df1 <- 1
  df2 <- length(x) - 2
  F <- (ssb / df1) / (ssw / df2)
  p <- pf(F, df1, df2, lower.tail = FALSE)
  data.frame(
    Riset = nama_riset,
    Uji = "Levene (median)",
    F = F,
    df1 = df1,
    df2 = df2,
    `p-value` = p,
    Keputusan = if (p < alpha) "Varians tidak homogen" else "Varians homogen",
    check.names = FALSE
  )
}

uji_fmax <- function(g1, g2, nama_riset) {
  v <- c(var(g1), var(g2))
  fmax <- max(v) / min(v)
  data.frame(
    Riset = nama_riset,
    Uji = "Hartley Fmax",
    Fmax = fmax,
    Keputusan = if (fmax > 10) "Indikasi tidak homogen" else "Indikasi homogen",
    check.names = FALSE
  )
}

tabel_varians <- function(g1, g2, nama_riset) {
  data.frame(
    Riset = nama_riset,
    Kelompok = c("Kelompok 1", "Kelompok 2"),
    Varians = c(var(g1), var(g2)),
    SD = c(sd(g1), sd(g2)),
    check.names = FALSE
  )
}

plot_hist <- function(g1, g2, judul) {
  d <- data.frame(
    nilai    = c(g1, g2),
    kelompok = rep(c("Kelompok 1", "Kelompok 2"), c(length(g1), length(g2)))
  )
  ggplot(d, aes(x = nilai, fill = kelompok)) +
    geom_histogram(bins = 20, alpha = 0.5, position = "identity", color = "white") +
    labs(title = judul, x = "Nilai", y = "Jumlah", fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
}

plot_box <- function(g1, g2, judul) {
  d <- data.frame(
    nilai    = c(g1, g2),
    kelompok = rep(c("Kelompok 1", "Kelompok 2"), c(length(g1), length(g2)))
  )
  ggplot(d, aes(x = kelompok, y = nilai, fill = kelompok)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
    labs(title = judul, x = NULL, y = "Nilai") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}

ui <- fluidPage(
  titlePanel("Simulasi Uji t Independen: Dua Riset"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Pengaturan umum"),
      numericInput("n", "Ukuran sampel per kelompok (N)", value = 30, min = 3, max = 1000, step = 1),
      hr(),
      h4("Riset A"),
      h5("Kelompok 1"),
      numericInput("a_m1", "Mean", value = 10, step = 0.1),
      numericInput("a_sd1", "SD", value = 1, min = 0.01, step = 0.1),
      h5("Kelompok 2"),
      numericInput("a_m2", "Mean", value = 11, step = 0.1),
      numericInput("a_sd2", "SD", value = 1, min = 0.01, step = 0.1),
      hr(),
      h4("Riset B"),
      h5("Kelompok 1"),
      numericInput("b_m1", "Mean", value = 10, step = 0.1),
      numericInput("b_sd1", "SD", value = 1, min = 0.01, step = 0.1),
      h5("Kelompok 2"),
      numericInput("b_m2", "Mean", value = 11, step = 0.1),
      numericInput("b_sd2", "SD", value = 1, min = 0.01, step = 0.1),
      hr(),
      checkboxInput("var_equal", "Asumsikan varians sama (Student t)", value = FALSE),
      numericInput("alpha", "Taraf signifikansi (alpha)", value = 0.05, min = 0.001, max = 0.5, step = 0.005),
      numericInput("seed", "Seed dasar", value = 123, step = 1),
      actionButton("simulasi", "Jalankan simulasi", icon = icon("refresh"), class = "btn-primary"),
      helpText("Ketik angka, lalu tekan tombol ini. Nilai akan dikunci pada saat tombol ditekan.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Perbandingan kedua riset",
          h3("Parameter aktif"),
          h4("Pengaturan umum"),
          tableOutput("param_umum"),
          h4("Riset A"),
          tableOutput("param_A"),
          h4("Riset B"),
          tableOutput("param_B"),
          hr(),
          h3("Perbandingan hasil uji t"),
          tableOutput("ringkas")
        ),
        tabPanel("Hasil analisis",
          h3("Riset A"),
          h4("Statistik deskriptif"),
          tableOutput("a_desk"),
          h4("Hasil uji t"),
          tableOutput("a_uji"),
          hr(),
          h3("Riset B"),
          h4("Statistik deskriptif"),
          tableOutput("b_desk"),
          h4("Hasil uji t"),
          tableOutput("b_uji")
        ),
        tabPanel("Plot",
          h3("Riset A"),
          plotOutput("a_hist", width = "650px", height = "320px"),
          plotOutput("a_box", width = "450px", height = "320px"),
          hr(),
          h3("Riset B"),
          plotOutput("b_hist", width = "650px", height = "320px"),
          plotOutput("b_box", width = "450px", height = "320px")
        ),
        tabPanel("Uji asumsi",
          h3("Uji normalitas (Shapiro-Wilk)"),
          tableOutput("norm"),
          hr(),
          h3("Uji homogenitas varians (Levene)"),
          tableOutput("levene"),
          hr(),
          h3("Hartley Fmax"),
          tableOutput("fmax"),
          hr(),
          h3("Varians per kelompok"),
          tableOutput("varian")
        ),
        tabPanel("Data",
          h3("Data hasil bangkitan"),
          p("Setiap klik tombol Jalankan simulasi akan membangkitkan data baru."),
          h4("Format lebar"),
          tableOutput("data_lebar"),
          h4("Format panjang"),
          tableOutput("data_panjang"),
          hr(),
          h4("Unduh"),
          downloadButton("unduh_data", "Unduh data (CSV)"),
          downloadButton("unduh_ringkas", "Unduh ringkasan hasil (CSV)")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  p <- reactiveValues(
    n = 30,
    a_m1 = 10, a_sd1 = 1, a_m2 = 11, a_sd2 = 1,
    b_m1 = 10, b_sd1 = 1, b_m2 = 11, b_sd2 = 1,
    var_equal = FALSE, alpha = 0.05, seed = 123, run = 0
  )

  observeEvent(input$simulasi, {
    p$n <- input$n
    p$a_m1 <- input$a_m1
    p$a_sd1 <- input$a_sd1
    p$a_m2 <- input$a_m2
    p$a_sd2 <- input$a_sd2
    p$b_m1 <- input$b_m1
    p$b_sd1 <- input$b_sd1
    p$b_m2 <- input$b_m2
    p$b_sd2 <- input$b_sd2
    p$var_equal <- input$var_equal
    p$alpha <- input$alpha
    p$seed <- input$seed
    p$run <- p$run + 1
  })

  sim <- reactive({
    validate(
      need(!is.na(p$n) && p$n >= 3, "N minimal 3."),
      need(!is.na(p$a_sd1) && !is.na(p$a_sd2) && !is.na(p$b_sd1) && !is.na(p$b_sd2) &&
             p$a_sd1 > 0 && p$a_sd2 > 0 && p$b_sd1 > 0 && p$b_sd2 > 0,
           "SD harus lebih besar dari 0."),
      need(!is.na(p$alpha) && p$alpha > 0 && p$alpha < 1, "Alpha harus di antara 0 dan 1.")
    )
    set.seed(p$seed + p$run)
    n <- as.integer(p$n)
    a1 <- rnorm(n, mean = p$a_m1, sd = p$a_sd1)
    a2 <- rnorm(n, mean = p$a_m2, sd = p$a_sd2)
    b1 <- rnorm(n, mean = p$b_m1, sd = p$b_sd1)
    b2 <- rnorm(n, mean = p$b_m2, sd = p$b_sd2)
    list(A = list(g1 = a1, g2 = a2), B = list(g1 = b1, g2 = b2))
  })

  data_lebar <- reactive({
    s <- sim()
    data.frame(
      id   = seq_len(as.integer(p$n)),
      A_K1 = s$A$g1,
      A_K2 = s$A$g2,
      B_K1 = s$B$g1,
      B_K2 = s$B$g2,
      check.names = FALSE
    )
  })

  data_panjang <- reactive({
    s <- sim()
    n <- as.integer(p$n)
    rbind(
      data.frame(Riset = "Riset A", Kelompok = "Kelompok 1", id = seq_len(n), nilai = s$A$g1),
      data.frame(Riset = "Riset A", Kelompok = "Kelompok 2", id = seq_len(n), nilai = s$A$g2),
      data.frame(Riset = "Riset B", Kelompok = "Kelompok 1", id = seq_len(n), nilai = s$B$g1),
      data.frame(Riset = "Riset B", Kelompok = "Kelompok 2", id = seq_len(n), nilai = s$B$g2)
    )
  })

  hasil_A <- reactive({
    s <- sim()$A
    hitung_uji(s$g1, s$g2, p$var_equal, p$alpha)
  })

  hasil_B <- reactive({
    s <- sim()$B
    hitung_uji(s$g1, s$g2, p$var_equal, p$alpha)
  })

  output$param_umum <- renderTable({
    data.frame(
      Parameter = c("N per kelompok", "Alpha", "Uji varians", "Seed"),
      Nilai = c(p$n, p$alpha,
                if (isTRUE(p$var_equal)) "Student" else "Welch",
                p$seed)
    )
  })

  output$param_A <- renderTable({
    data.frame(
      Kelompok = c("Kelompok 1", "Kelompok 2"),
      Mean = c(p$a_m1, p$a_m2),
      SD = c(p$a_sd1, p$a_sd2)
    )
  })

  output$param_B <- renderTable({
    data.frame(
      Kelompok = c("Kelompok 1", "Kelompok 2"),
      Mean = c(p$b_m1, p$b_m2),
      SD = c(p$b_sd1, p$b_sd2)
    )
  })

  output$ringkas <- renderTable({
    rbind(
      cbind(Riset = "Riset A", hasil_A()$ringkas),
      cbind(Riset = "Riset B", hasil_B()$ringkas)
    )
  }, digits = 4)

  output$a_desk <- renderTable({
    s <- sim()$A
    rbind(
      desk(s$g1, "Kelompok 1", p$alpha, p$a_m1, p$a_sd1),
      desk(s$g2, "Kelompok 2", p$alpha, p$a_m2, p$a_sd2)
    )
  }, digits = 3)

  output$b_desk <- renderTable({
    s <- sim()$B
    rbind(
      desk(s$g1, "Kelompok 1", p$alpha, p$b_m1, p$b_sd1),
      desk(s$g2, "Kelompok 2", p$alpha, p$b_m2, p$b_sd2)
    )
  }, digits = 3)

  output$a_uji <- renderTable({
    hasil_A()$tabel
  })

  output$b_uji <- renderTable({
    hasil_B()$tabel
  })

  output$norm <- renderTable({
    s <- sim()
    rbind(
      uji_normalitas(s$A$g1, s$A$g2, "Riset A", p$alpha),
      uji_normalitas(s$B$g1, s$B$g2, "Riset B", p$alpha)
    )
  }, digits = 4)

  output$levene <- renderTable({
    s <- sim()
    rbind(
      uji_levene(s$A$g1, s$A$g2, "Riset A", p$alpha),
      uji_levene(s$B$g1, s$B$g2, "Riset B", p$alpha)
    )
  }, digits = 4)

  output$fmax <- renderTable({
    s <- sim()
    rbind(
      uji_fmax(s$A$g1, s$A$g2, "Riset A"),
      uji_fmax(s$B$g1, s$B$g2, "Riset B")
    )
  }, digits = 4)

  output$varian <- renderTable({
    s <- sim()
    rbind(
      tabel_varians(s$A$g1, s$A$g2, "Riset A"),
      tabel_varians(s$B$g1, s$B$g2, "Riset B")
    )
  }, digits = 4)

  output$data_lebar <- renderTable({
    data_lebar()
  }, digits = 4)

  output$data_panjang <- renderTable({
    data_panjang()
  }, digits = 4)

  output$a_hist <- renderPlot({
    s <- sim()$A
    plot_hist(s$g1, s$g2, "Riset A: histogram")
  })

  output$a_box <- renderPlot({
    s <- sim()$A
    plot_box(s$g1, s$g2, "Riset A: boxplot")
  })

  output$b_hist <- renderPlot({
    s <- sim()$B
    plot_hist(s$g1, s$g2, "Riset B: histogram")
  })

  output$b_box <- renderPlot({
    s <- sim()$B
    plot_box(s$g1, s$g2, "Riset B: boxplot")
  })

  output$unduh_data <- downloadHandler(
    filename = function() paste0("data_simulasi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      write.csv(data_panjang(), file, row.names = FALSE)
    }
  )

  output$unduh_ringkas <- downloadHandler(
    filename = function() paste0("ringkasan_hasil_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      s <- sim()
      d <- rbind(
        cbind(Riset = "Riset A", desk(s$A$g1, "Kelompok 1", p$alpha, p$a_m1, p$a_sd1)),
        cbind(Riset = "Riset A", desk(s$A$g2, "Kelompok 2", p$alpha, p$a_m2, p$a_sd2)),
        cbind(Riset = "Riset B", desk(s$B$g1, "Kelompok 1", p$alpha, p$b_m1, p$b_sd1)),
        cbind(Riset = "Riset B", desk(s$B$g2, "Kelompok 2", p$alpha, p$b_m2, p$b_sd2))
      )
      u <- rbind(
        cbind(Riset = "Riset A", hasil_A()$ringkas),
        cbind(Riset = "Riset B", hasil_B()$ringkas)
      )
      write.csv(d, file, row.names = FALSE)
      write.table(u, file, sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE)
    }
  )
}

shinyApp(ui, server)
