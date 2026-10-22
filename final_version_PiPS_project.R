library(shiny)
library(ggplot2)

# Helper: both ANOVA designs use F tests
is_f <- function(test) test %in% c('anova', 'anova2')

# Power function
# 'design' carries the extra design information for the ANOVA tests:
#   one-way ANOVA ('anova'):   design = k (number of groups)
#   two-way ANOVA ('anova2'):  design = list(a = levels of A, b = levels of B, effect = 'A' | 'B' | 'AB')
# For ANOVA, n is the sample size per group (one-way) or per cell (two-way); d is Cohen's f.

calc_power <- function(n, d, alpha, test, alternative, design = 3) {
  if (is_f(test)) {
    if (test == 'anova') {
      cells <- design
      df1 <- cells - 1
    } else {
      cells <- design$a * design$b
      df1 <- switch(design$effect,
                    A  = design$a - 1,
                    B  = design$b - 1,
                    AB = (design$a - 1) * (design$b - 1))
    }
    N <- cells * n
    df2 <- N - cells
    lambda <- d^2 * N
    cut <- qf(1 - alpha, df1, df2)
    power <- pf(cut, df1, df2, ncp = lambda, lower.tail = FALSE)
    return(list(power = power, df = df2, df1 = df1, delta = lambda,
                critical = cut, total = N))
  }
  # One and two-sample t-tests
  df <- if (test == 'one') n - 1 else 2 * n - 2
  delta <- if (test == 'one') d * sqrt(n) else d * sqrt(n / 2)
  if (alternative == 'two') {
    cut <- qt(1 - alpha / 2, df)
    power <- pt(-cut, df, ncp = delta) + pt(cut, df, ncp = delta, lower.tail = FALSE)
    critical <- c(-cut, cut)
  } else if (alternative == 'greater') {
    cut <- qt(1 - alpha, df)
    power <- pt(cut, df, ncp = delta, lower.tail = FALSE)
    critical <- cut
  } else {
    cut <- qt(alpha, df)
    power <- pt(cut, df, ncp = delta)
    critical <- cut
  }
  list(power = power, df = df, delta = delta, critical = critical,
       total = if (test == 'one') n else 2 * n)
}

# Density helpers and rejection region that work for both t and F tests
dens0 <- function(x, r, test) if (is_f(test)) df(x, r$df1, r$df) else dt(x, r$df)
# For negative effects, dt() with a negative ncp emits precision warnings; the mirror identity
# dt(x, df, -δ) = dt(-x, df, δ) gives the same density without them.
dens1 <- function(x, r, test) {
  if (is_f(test)) return(df(x, r$df1, r$df, ncp = r$delta))
  s <- if (r$delta < 0) -1 else 1
  dt(s * x, r$df, ncp = abs(r$delta))
}
reject_region <- function(x, r, alternative, test) {
  if (is_f(test)) x >= r$critical
  else if (alternative == 'two') abs(x) >= max(abs(r$critical))
  else if (alternative == 'greater') x >= r$critical
  else x <= r$critical
}
stat_label <- function(test) if (is_f(test)) 'F-statistic' else 't-statistic'

# Upper y-limit for the density plots. An F density with 1 numerator df is infinite at 0,
# so for F tests the very left edge of the curve is ignored when choosing the axis height.
y_top <- function(dat, r, test) {
  keep <- if (is_f(test)) dat$x >= 0.1 * r$critical else rep(TRUE, nrow(dat))
  max(dat$h0[keep], dat$h1[keep], na.rm = TRUE) * 1.12
}



# Help texts

help_text <- list(
  test = c('Which test?', 'One-sample t-test: one mean vs a fixed value. Two-sample t-test: the means of two independent groups. One-way ANOVA: the means of three or more groups. Two-way ANOVA: two factors and their interaction. All assume independent, normally distributed observations; the group comparisons also assume equal variances and, in this app, equal group sizes.'),
  hypothesis = c('Hypotheses (H₀ and H₁)', 'H₀ says there is no effect; H₁ says there is one. Power is calculated for the specific effect you enter. Not rejecting H₀ does not prove it is true.'),
  pvalue = c('Significance and p-values', 'The p-value is the probability of a result at least this extreme if H₀ is true. If p ≤ α the result is significant. It is not the probability that H₀ is true, and it does not measure the size of the effect.'),
  alternative = c('Alternative hypothesis', 'A two-sided test rejects H₀ for a difference in either direction; a one-sided test (Greater or Less) only in the specified direction. Choose before seeing the data.'),
  n = c('Sample size (n)', 'Participants per group (per cell in the two-way ANOVA). For the one-sample test, the total number of observations.'),
  d = c("Effect size", 'Cohen’s d (t-tests): the difference between the means divided by the standard deviation. Cohen’s f (ANOVA): the standard deviation of the group means divided by the within-group standard deviation. Conventional benchmarks: d = 0.2 / 0.5 / 0.8 and f = 0.10 / 0.25 / 0.40 for small / medium / large.'),
  alpha = c('Significance level (α)', 'The probability of rejecting H₀ when it is true (a Type I error). Usually set at .05.'),
  target = c('Target power', 'The power you want to reach. 80% is a common choice.'),
  power = c('Statistical power', 'The probability of rejecting H₀ when the assumed effect is real (1 − β).'),
  total = c('Total participants', 'All participants across groups (or cells).'),
  required = c('Required n', 'The smallest sample size n (per group or cell) that reaches your target power.'),
  df = c('Degrees of freedom', 'Determine the exact shape of the t- or F-distribution. They depend on the sample size and the number of groups.'),
  ncp = c('Noncentrality parameter', 'How far the H₁ distribution is shifted away from H₀. It grows with both effect size and sample size.'),
  critical = c('Critical value', 'The cut-off for significance. Results beyond it reject H₀.'),
  null = c('Blue curve: H₀', 'The sampling distribution of the test statistic when H₀ is true.'),
  alt = c('Pink curve: H₁', 'The sampling distribution of the test statistic when the assumed effect is real.'),
  region = c('Purple area: power', 'The area under the H₁ curve beyond the critical value(s). This area is the power.'),
  boundary = c('Dashed line: critical value', 'Results beyond this line are significant.'),
  curve = c('Power curve', 'Power for each sample size, with all other settings unchanged.'),
  targetline = c('Target line', 'Your target power. Where the curve crosses it is the required n.'),
  current = c('Your design', 'Your current sample size and its power.'),
  dist_note = c('Explore the distributions', 'Blue: no effect (H₀). Pink: your assumed effect (H₁). The purple area beyond the dashed line is the power. Click any part of the chart for details.'),
  errors_note = c('Understand α, β and power', 'Power (1 − β): the probability of rejecting H₀ when H₁ is true. Type II error (β): the probability of not rejecting H₀ when H₁ is true. Type I error (α): the probability of rejecting H₀ when H₀ is true. 1 − α: the probability of not rejecting H₀ when H₀ is true.'),
  curve_note = c('Power playground', 'Shows how power grows with sample size. The dashed line is your target; the dot is your design.'),
  k = c('Number of groups (k)', 'How many groups the one-way ANOVA compares.'),
  levels_a = c('Levels of A', 'How many categories factor A has.'),
  levels_b = c('Levels of B', 'How many categories factor B has.'),
  effect = c('Effect to test', 'Main A or Main B: whether the means differ between the levels of that factor, averaged over the other factor. A × B: whether the effect of A differs across the levels of B. Enter Cohen’s f for the effect you choose.')
)


info_button <- function(text, key) {
  div(class = 'labelrow', tags$label(text),
      actionButton(paste0('help_', key), 'ⓘ', class = 'info', title = paste('Learn about', text)))
}

# Card heading with an ⓘ button; the explanation lives in help_text[[key]] instead of on the page.
card_title <- function(title, key) {
  div(class = 'card-title', h3(title),
      actionButton(paste0('help_', key), 'ⓘ', class = 'info', title = paste('About:', title)))
}

# A slider with a number box next to it. The box has the id '<id>_num'; the server keeps both in sync,
# and all calculations keep reading the slider value (input[[id]]).
slider_limits <- list(n = c(3, 200, 1), k = c(3, 10, 1), levels_a = c(2, 6, 1), levels_b = c(2, 6, 1),
                      d = c(-2, 2, .01), alpha = c(.001, .10, .001), target = c(.7, .95, .01))
slider_box <- function(id, min, max, value, step) {
  div(class = 'param-row',
      div(class = 'param-slider', sliderInput(id, NULL, min, max, value, step = step, ticks = FALSE, width = '100%')),
      div(class = 'param-box', tagAppendAttributes(
        numericInput(paste0(id, '_num'), NULL, value, min = min, max = max, step = step, width = '100%'),
        `aria-label` = paste('Type a value for', id), .cssSelector = 'input')))
}

# UI

ui <- fluidPage(
  tags$head(
    tags$meta(name = 'viewport', content = 'width=device-width, initial-scale=1'),
    # Charts on the Explore page stretch to the available height. Shiny only re-renders a plot when
    # the window resizes, so tell it to re-measure whenever a chart box changes size.
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        if (!window.ResizeObserver) return;
        var timer;
        var ro = new ResizeObserver(function() {
          clearTimeout(timer);
          timer = setTimeout(function() { $(window).trigger('resize'); }, 120);
        });
        document.querySelectorAll('.plot-card .shiny-plot-output').forEach(function(el) { ro.observe(el); });
      });")),
    tags$style(HTML(" 
      @import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
      :root{--bg:#090f21;--panel:#121b32;--panel2:#182440;--edge:#2a3753;--text:#edf2ff;--muted:#a4b1ce;--violet:#9476ff;--pink:#ff729d;--blue:#58b8ff}
      body{background:var(--bg);color:var(--text);font-family:'DM Sans',system-ui,sans-serif;font-size:15px}
      h1,h2,h3,h4,.brand{font-family:'Space Grotesk',system-ui,sans-serif;letter-spacing:-.035em}
      .container-fluid{max-width:1450px;padding:20px 28px}
      .navbar{display:none}
      .topbar{display:flex;align-items:center;justify-content:space-between;margin:0 0 24px;gap:14px}
      .brand{font-weight:700;font-size:27px;color:#fff}.brand span{color:var(--violet)}
      .subbrand{font-size:13px;color:var(--muted)}
      .steps{display:flex;gap:8px;align-items:center;flex-wrap:wrap;color:var(--muted);font-size:13px}
      .step{background:#192641;border:1px solid var(--edge);padding:9px 13px;border-radius:28px}
      .step.active{background:#352960;color:#e8deff;border-color:#7865dc}
      .cardx{background:var(--panel);border:1px solid var(--edge);border-radius:19px;padding:22px;margin-bottom:16px}
      .cardx h2,.cardx h3{margin:0 0 12px}.cardx h3{font-size:20px}
      .eyebrow{color:#bcaeff;text-transform:uppercase;letter-spacing:.14em;font-weight:700;font-size:12px}
      .hero{padding:50px 42px;background:radial-gradient(circle at 85% 10%,#33285e 0,transparent 52%),linear-gradient(125deg,#182347,#10192e);border:1px solid #40416f;border-radius:26px;margin:12px 0 22px}
      .hero h1{font-size:clamp(38px,5vw,68px);font-weight:700;margin:12px 0}.hero h1 span{color:#aa88ff}
      .hero p{color:#cad3ed;max-width:680px;font-size:18px;line-height:1.65}
      .btn-primary,.btn-start{background:linear-gradient(110deg,#6657e9,#a05cf4);border:0;border-radius:11px;padding:12px 22px;color:#fff;font-weight:700;transition:transform .2s,filter .2s}
      .btn-primary:hover,.btn-start:hover{color:#fff;filter:brightness(1.12);transform:translateY(-2px)}
      .btn-secondary{border:1px solid #65658e;background:#222d4d;color:white;border-radius:11px;padding:10px 17px}
      .labelrow{display:flex;align-items:center;gap:6px;margin-top:17px;margin-bottom:6px}
      .labelrow label{margin:0;font-size:13px;font-weight:700;color:#dbe3fc}
      .info{border:0;background:none;color:#aa96ff;font-size:16px;font-weight:700;padding:2px 6px;box-shadow:none}
      .info:hover,.info:focus{background:#393158;color:#fff;border-radius:9px}
      .form-control,.selectize-input,.selectize-dropdown{background:#202c49!important;color:#fff!important;border:1px solid #455271!important;border-radius:9px!important}
      .radio label{color:#dce5f9;margin-bottom:7px}
      .irs--shiny .irs-bar,.irs--shiny .irs-single{background:#946bff;border-color:#946bff}
      .irs--shiny .irs-handle{border-color:#9675ff;background:#d9cbff}
      .irs--shiny .irs-min,.irs--shiny .irs-max{background:transparent;color:#98a9cb}
      .irs--shiny .irs-grid-text,.irs--shiny .irs-from,.irs--shiny .irs-to{color:#d2ddfb}
      .metric{font:700 38px 'Space Grotesk',sans-serif;color:#b19aff}.metric-sub{font-size:13px;color:var(--muted);line-height:1.5}
      .helper{font-size:13px;color:var(--muted);line-height:1.7}
      .chart-note{font-size:12px;color:#abb8d1}
      .chart-key{display:flex;flex-wrap:wrap;gap:8px;margin:8px 0}
      .chart-key .btn{background:#202c4b;color:#e3eaff;border:1px solid #405070;border-radius:20px;font-size:12px;padding:7px 12px;white-space:normal}
      .chart-key .btn:hover{border-color:#a48cff;background:#343154}
      .plotly,.shiny-plot-output{max-width:100%}
      .modal-content{background:#192540;color:#edf2ff;border:1px solid #7560be;border-radius:20px;box-shadow:0 24px 80px #000a;animation:modalIn .23s ease both;overflow:hidden}
      .modal-header,.modal-footer{border-color:#354260}.modal-title{font-family:'Space Grotesk';font-weight:700}.modal-body{font-size:15px;line-height:1.8;color:#d7e0f5}
      .modal-footer .btn{background:#6858db;color:white;border:0;border-radius:9px;padding:9px 18px}
      .modal-backdrop.in{opacity:.7}@keyframes modalIn{from{transform:translateY(12px) scale(.975);opacity:0}to{transform:translateY(0) scale(1);opacity:1}}
      .modebadge{display:inline-block;background:#392e64;color:#e0d5ff;padding:7px 12px;border-radius:12px;font-weight:700} .mini{background:#1a2847;border:1px solid #45537a;border-radius:15px;padding:18px;margin:10px 0;min-height:125px}.intro-graphic{height:110px;display:flex;align-items:flex-end;gap:8px;margin:22px 0}.intro-graphic span{display:block;width:12%;background:linear-gradient(#bf8bff,#565fe8);border-radius:8px 8px 0 0;animation:grow .65s ease both}@keyframes grow{from{height:4px}} .mode-note{color:#c8d3f2;font-size:13px}.learn-tab .btn{margin-right:8px;margin-bottom:8px} .text-warning{color:#ffba86!important} .footnote{color:#8b9bbf;font-size:12px;margin:20px 4px 10px}
      @media(max-width:767px){.container-fluid{padding:15px}.hero{padding:25px}.hero p{font-size:16px}.metric{font-size:30px}.steps{justify-content:flex-start}}

      /* PowerLab 3.0 — consistent visual language and spacious type scale */
      body{font-size:17px;line-height:1.72;-webkit-font-smoothing:antialiased}
      p{line-height:1.8;margin:0 0 18px}
      .container-fluid{max-width:1450px;padding:32px 34px}
      .cardx{padding:clamp(24px,3vw,38px);margin-bottom:24px;border-radius:24px;
        box-shadow:0 12px 34px rgba(0,0,0,.10)}
      .cardx h2{font-size:clamp(25px,3vw,34px);line-height:1.22;margin-bottom:19px}
      .cardx h3{font-size:23px;line-height:1.3;margin-bottom:22px}
      .cardx h4{font-size:20px;line-height:1.35}
      .hero{padding:clamp(34px,6vw,84px);min-height:520px;display:flex;flex-direction:column;
        justify-content:center;gap:12px;margin-top:30px;border-radius:30px;
        background:radial-gradient(circle at 82% 22%,rgba(151,112,255,.23),transparent 43%),
          radial-gradient(circle at 10% 90%,rgba(88,184,255,.10),transparent 36%),
          linear-gradient(130deg,#192849,#10172d)}
      .hero h1{font-size:clamp(44px,6.3vw,84px);line-height:1.08;max-width:870px;
        letter-spacing:-.06em;margin:16px 0 10px}
      .hero p{font-size:clamp(18px,2vw,22px);line-height:1.78;max-width:760px;margin-bottom:12px}
      .hero .eyebrow{font-size:13px;letter-spacing:.2em}
      .section-copy{font-size:18px;line-height:1.88;color:#cbd5ed;max-width:850px;margin:0 0 28px}
      .stepper{display:flex;flex-wrap:wrap;gap:11px;margin:12px 0 27px}
      .step-pill{display:inline-flex;align-items:center;padding:11px 16px;
        border:1px solid #334361;background:#141f37;border-radius:100px;
        color:#aab9d6;font-size:13px;line-height:1.4;font-weight:600}
      .step-pill.active{color:#f4efff;border-color:#9476ff;background:#332958;
        box-shadow:0 0 0 3px rgba(148,118,255,.10)}
      .learn-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px;
        margin:27px 0 34px}
      .learn-card{position:relative;display:flex;flex-direction:column;align-items:flex-start;
        min-height:267px;padding:28px 26px;gap:8px;border-radius:22px;
        background:linear-gradient(150deg,#1e2a49,#17223c);
        border:1px solid #374764;box-shadow:0 10px 22px rgba(0,0,0,.12);
        transition:border-color .22s,transform .22s,box-shadow .22s}
      .learn-card:hover{border-color:#8774d7;transform:translateY(-3px);
        box-shadow:0 16px 28px rgba(0,0,0,.18)}
      .learn-card:before{content:'✦';display:grid;place-items:center;width:42px;height:42px;
        margin-bottom:7px;border-radius:13px;color:#d9c8ff;
        background:linear-gradient(135deg,#49367b,#2a3568);font-size:21px}
      .learn-card:nth-child(2):before{content:'H₀';font-size:18px}
      .learn-card:nth-child(3):before{content:'p';font-style:italic}
      .learn-card:nth-child(4):before{content:'α'}
      .learn-card:nth-child(5):before{content:'d'}
      .learn-card:nth-child(6):before{content:'1−β';font-size:14px}
      .learn-card h4{color:#f4f1ff;margin:0 0 5px}
      .learn-card p{font-size:15.5px;color:#c3d0e9;line-height:1.78;margin:0 0 16px}
      .learn-card .btn{margin-top:auto;font-size:14px;padding:9px 15px;border-radius:10px}
      .choice-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:17px;
        margin:24px 0 30px}
      .choice-card.btn{width:100%;min-height:150px;padding:24px 24px 24px 26px;
        display:flex;flex-direction:column;justify-content:center;align-items:flex-start;
        position:relative;text-align:left;white-space:normal;
        border:1px solid #394965;border-radius:21px;color:#f1efff;
        background:linear-gradient(145deg,#1e2a48,#162139);
        box-shadow:0 10px 25px rgba(0,0,0,.10);
        transition:transform .2s,border-color .2s,background .2s,box-shadow .2s}
      .choice-card.btn:hover,.choice-card.btn:focus-visible{
        color:#fff;transform:translateY(-3px);border-color:#a28aff;
        background:linear-gradient(145deg,#29365a,#231f45);
        box-shadow:0 12px 32px rgba(112,79,205,.18);outline:none}
      .choice-card.btn.selected{border:2px solid #aa96ff;
        background:linear-gradient(145deg,#352b61,#202947);
        box-shadow:0 0 0 3px rgba(148,118,255,.12)}
      .choice-name{display:block;font-family:'Space Grotesk',sans-serif;font-weight:700;
        font-size:18px;line-height:1.4;margin-bottom:12px;white-space:normal}
      .choice-description{display:block;color:#c1cce5;font-size:15px;line-height:1.65;
        font-weight:400;white-space:normal}
      .choice-card.test-card.btn{min-height:118px;padding:21px 22px}
      .choice-card.test-card .choice-name{font-size:16px}
      .choice-card.test-card .choice-description{font-size:14px}
      .explain-callout{padding:25px 28px;border:1px solid #64549e;border-left:4px solid #ab8bff;
        border-radius:17px;background:linear-gradient(110deg,#29264c,#202a47);
        font-size:17px;line-height:1.85;color:#e3dffb;margin:25px 0 30px}
      .explain-callout strong{color:#fff}
      .nav-actions{display:flex;align-items:center;flex-wrap:wrap;gap:14px;margin-top:28px}
      .btn-start,.btn-primary{font-size:16px;line-height:1.4;padding:14px 23px}
      .btn-secondary{font-size:15px;line-height:1.4;padding:12px 19px}
      .labelrow{margin-top:25px;margin-bottom:12px}
      .labelrow label{font-size:15px}
      .helper,.chart-note,.plot-help{font-size:15px;line-height:1.76;color:#b6c4df}
      .metric-sub,.footnote{font-size:14px;line-height:1.75}
      .modal-body{font-size:17px;line-height:1.85;padding:24px 29px}
      .modal-header{padding:23px 29px}
      .modal-title{font-size:23px}
      .modal-footer{padding:16px 29px}
      .radio label{line-height:1.65}
      .chart-key{gap:10px;margin:14px 0 18px}
      .chart-key .btn{font-size:14px;line-height:1.45;padding:9px 14px}
      .shiny-input-container{margin-bottom:20px}
      .topbar{margin:0 0 24px}
      @media(max-width:1100px){.learn-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
      @media(max-width:767px){
        body{font-size:16px}.container-fluid{padding:18px 14px}
        .cardx{padding:22px 18px;margin-bottom:18px}
        .hero{padding:35px 23px;min-height:480px;margin-top:9px}
        .hero h1{font-size:clamp(39px,9vw,57px)}
        .hero p{font-size:17px}.section-copy{font-size:16px}
        .learn-grid,.choice-grid{grid-template-columns:1fr;gap:13px}
        .learn-card{min-height:0;padding:23px}
        .choice-card.btn{min-height:125px}
        .step-pill{font-size:12px;padding:9px 11px}
        .nav-actions .btn{width:100%;white-space:normal}
      }


      /* Compact controls in the main lab; onboarding cards remain large. */
      .test-choice-list{display:flex;flex-direction:column;gap:9px;margin:15px 0 18px}
      .test-choice-list .choice-card.test-card.btn{width:100%;min-height:48px;height:auto;
        padding:12px 15px;border-radius:12px;display:flex;flex-direction:row;
        align-items:center;justify-content:space-between;gap:12px;
        text-align:left;box-shadow:none;transform:none}
      .test-choice-list .choice-card.test-card.btn:hover,
      .test-choice-list .choice-card.test-card.btn:focus-visible{transform:none;box-shadow:0 0 0 2px rgba(148,118,255,.15)}
      .test-choice-list .choice-card.test-card .choice-name{font-family:'DM Sans',sans-serif;
        font-size:14px;font-weight:600;line-height:1.45;margin:0;white-space:normal}
      .test-choice-list .test-check{color:#d4c5ff;min-width:15px;text-align:right;font-weight:700}
      .test-modal-actions{display:flex;justify-content:flex-end;gap:12px;flex-wrap:wrap;width:100%}
      .test-modal-actions .btn{min-height:44px;min-width:150px;padding:11px 18px;
        font-family:'DM Sans',sans-serif;font-size:15px;font-weight:700;
        line-height:1.4;border-radius:11px;display:inline-flex;align-items:center;justify-content:center}
      .test-modal-actions .btn-default{background:#273556;color:#edf2ff;border:1px solid #63708d}
      .test-modal-actions .btn-default:hover{color:#fff;background:#364669;border-color:#a48cff}
      .test-modal-actions .btn-primary{background:linear-gradient(110deg,#6657e9,#a05cf4);border:0}
      @media(max-width:480px){.test-modal-actions .btn{flex:1 1 120px;min-width:0}}


      /* UX MERGE: compact test selector with separate context-help buttons. */
      .test-select-row{display:flex;align-items:stretch;gap:9px;width:100%}
      .test-select-row .choice-card.test-card.btn{flex:1;min-width:0;width:auto}
      .test-select-row .test-info-button.btn{flex:0 0 48px;min-width:48px;
        border:1px solid #405070;border-radius:12px;background:#202c49;
        font-size:19px;line-height:1;display:flex;align-items:center;justify-content:center}
      .test-select-row .test-info-button.btn:hover,.test-select-row .test-info-button.btn:focus-visible{
        border-color:#a48cff;background:#393158;color:#fff}


      /* FINAL UX: readable onboarding, large preview input and polished two-level switch. */
      .cardx > .section-copy{max-width:none;width:100%;padding-right:0}
      .hero p{max-width:900px}
      .intro-preview-slider{width:100%;max-width:none;margin:26px 0 24px;padding:22px 25px 27px;
        background:linear-gradient(125deg,#202d4f,#192540);border:1px solid #46577d;border-radius:20px}
      .intro-preview-slider .shiny-input-container{width:100%!important;max-width:none;margin-bottom:0}
      .intro-preview-slider .control-label{font-size:18px;font-weight:700;margin-bottom:25px;color:#f1efff}
      .intro-preview-slider .irs{height:66px}
      .intro-preview-slider .irs--shiny .irs-line{height:9px;top:32px;border-radius:9px;background:#354569;border:0}
      .intro-preview-slider .irs--shiny .irs-bar{height:9px;top:32px;border-radius:9px}
      .intro-preview-slider .irs--shiny .irs-handle{width:19px;height:19px;top:27px;border:3px solid #d8caff;
        background:#8e75ff;box-shadow:0 0 0 3px rgba(149,118,255,.15);cursor:grab}
      .intro-preview-slider .irs--shiny .irs-single{font-size:16px;top:-4px;padding:4px 9px;border-radius:8px}
      /* Grid already labels endpoints (5 and 200); hide duplicate ionRangeSlider endpoint labels. */
      .intro-preview-slider .irs--shiny .irs-min,.intro-preview-slider .irs--shiny .irs-max{display:none!important}
      .intro-preview-slider .irs--shiny .irs-from,.intro-preview-slider .irs--shiny .irs-to{font-size:16px}
      .lab-hero{position:relative;overflow:hidden;padding:28px 32px;margin-bottom:24px;
        border:1px solid #524a85;border-radius:24px;
        background:radial-gradient(circle at 91% 2%,rgba(165,111,255,.24),transparent 43%),
          linear-gradient(120deg,#202d52,#151d36 78%);box-shadow:0 12px 36px #0002}
      .lab-hero-content{display:flex;justify-content:space-between;align-items:center;gap:26px;flex-wrap:wrap}
      .lab-hero .brand{font-size:31px;line-height:1.25;letter-spacing:-.05em;margin-bottom:7px}
      .lab-hero .lab-eyebrow{color:#c5aaff;text-transform:uppercase;letter-spacing:.14em;font-size:12px;font-weight:700;margin-bottom:7px}
      .lab-hero-actions{display:flex;align-items:center;gap:15px;flex-wrap:wrap}
      .level-switch{display:flex;gap:6px;background:#121b32;border:1px solid #495377;
        border-radius:14px;padding:5px}
      .level-switch .btn{font-size:15px;font-weight:700;border:1px solid transparent;
        background:transparent;color:#b9c6e2;border-radius:10px;padding:10px 20px;min-width:105px}
      .level-switch .btn:hover{background:#293359;color:#fff}
      .level-switch .btn.active{background:linear-gradient(120deg,#5f53cc,#8965eb);color:#fff;
        box-shadow:0 4px 13px #0003;border-color:#a78eff}
      .test-choice-list .test-check{display:none!important}
      @media(max-width:767px){.intro-preview-slider{padding:18px 15px 23px}
        .lab-hero{padding:23px 20px}.lab-hero-actions{width:100%;justify-content:space-between}
        .level-switch{width:100%}.level-switch .btn{flex:1;min-width:0}
        .lab-hero-actions>.btn{width:100%}}

      .warning-note{margin-top:16px;padding:13px 16px;border:1px solid #8a5a3a;border-left:4px solid #ffb776;
        border-radius:12px;background:#2c2333;color:#ffd9b8;font-size:14px;line-height:1.65}

      /* LAB PAGE: readable two-column layout. Parameters left; results and charts right. */
      .lab-page{font-size:16px;line-height:1.6}
      .lab-page .cardx{padding:20px 24px;margin:0;border-radius:20px}
      .lab-page .cardx h3{font-size:20px;line-height:1.3;margin:0 0 6px}
      .lab-page .lab-hero{padding:14px 24px;margin-bottom:16px;border-radius:20px}
      .lab-page .lab-hero .brand{font-size:27px;margin:0}
      .lab-page .lab-hero .lab-eyebrow{font-size:11px;margin-bottom:2px}
      .lab-page .level-switch .btn{padding:8px 18px;font-size:15px}
      .lab-page .labelrow{margin:16px 0 4px}.lab-page .labelrow label{font-size:15px}
      .lab-page .info{font-size:16px;padding:1px 6px}
      .lab-page .shiny-input-container{margin-bottom:0}
      .lab-page .radio-inline,.lab-page .radio label{font-size:15px;line-height:1.6}
      .lab-page .radio-inline{margin:0 18px 2px 0;padding-top:0}.lab-page .radio-inline+.radio-inline{margin-left:0}
      .lab-page .metric{font-size:34px;line-height:1.2;overflow-wrap:anywhere}
      .lab-page .chart-key{gap:8px;margin:6px 0 10px}
      .lab-page .chart-key .btn{font-size:13.5px;padding:5px 12px;line-height:1.4}
      .card-title{display:flex;align-items:center;gap:6px;justify-content:space-between}
      .card-title h3{margin:0!important;flex:1 1 auto}
      .btn-explain{background:#262f52;color:#dcd2ff;border:1px solid #5a4f99;border-radius:20px;
        font-size:14px;font-weight:600;padding:6px 14px}
      .btn-explain:hover,.btn-explain:focus{background:#393158;color:#fff;border-color:#a48cff}
      .error-caption{font-size:15px;color:#c8d3f2;margin-top:6px}
      .lab-left{display:flex;flex-direction:column;gap:16px;min-width:0}
      .lab-testcard .test-choice-list{margin:12px 0 4px;gap:10px}
      .lab-testcard .test-select-row{gap:8px}
      .lab-testcard .test-choice-list .choice-card.test-card.btn{min-height:54px;padding:12px 16px}
      .lab-testcard .test-choice-list .choice-card.test-card .choice-name{font-size:16px;line-height:1.35}
      .lab-testcard .test-select-row .test-info-button.btn{flex:0 0 50px;min-width:50px;font-size:20px}
      .param-pair{display:grid;grid-template-columns:1fr 1fr;gap:14px}
      .param-pair .param-box{flex-basis:60px}.param-pair .param-row{gap:8px}
      .lab-main{display:grid;grid-template-columns:minmax(330px,390px) minmax(0,1fr);gap:16px}
      .lab-params{min-width:0}
      .param-row{display:flex;align-items:flex-end;gap:12px}
      .param-slider{flex:1 1 auto;min-width:0}
      .param-box{flex:0 0 88px;padding-bottom:2px}
      .param-box .form-control{height:36px;padding:4px 8px;font-size:15px;text-align:center}
      .lab-right{display:flex;flex-direction:column;gap:16px;min-width:0}
      .metric-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}
      .lab-details{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px;font-size:15px;
        margin-top:6px;padding-top:6px;border-top:1px solid #2a3753}
      .lab-details .labelrow label{font-size:14px}
      /* Charts: in Advanced mode the two distribution charts are stacked (same x-axis, easy to compare),
         with the power curve beside them. In Beginner mode the power curve uses the full width. */
      .lab-plots{display:flex;gap:16px;align-items:flex-start}
      .dist-stack{flex:1.35 1 0;min-width:0;display:flex;flex-direction:column;gap:16px}
      .curve-cell{flex:1 1 0;min-width:0}
      .plot-card{min-width:0}
      @media(max-width:1099px){.lab-plots{flex-direction:column;align-items:stretch}}
      @media(max-width:899px){.lab-main{grid-template-columns:1fr}}

      @media(prefers-reduced-motion:reduce){*,*:before,*:after{animation:none!important;transition:none!important}}
    "))
  ),
   tabsetPanel(id = 'page', type = 'hidden',
       tabPanel(title='Welcome',value='intro',
         div(class='hero',div(class='eyebrow','WELCOME TO POWERLAB'),
             h1('Statistics you can ',span('see.')),
             p('Learn how researchers decide how many people to include in a study — and why some real effects are harder to detect than others.'),
             p('Explore statistical concepts through examples, interactive charts, and explanations designed for students.'),
             div(class='nav-actions',actionButton('intro_next1','Get started →',class='btn-start'),actionButton('intro_skip','Skip introduction',class='btn-secondary')))),
       tabPanel(title='The question',value='why',
         div(class='stepper',span(class='step-pill active','1 · Why power?'),span(class='step-pill','2 · Key concepts'),span(class='step-pill','3 · Tests'),span(class='step-pill','4 · Try it'),span(class='step-pill','5 · Your level')),
         div(class='cardx',h2('Why does statistical power matter?'),
             p(class='section-copy','Imagine you are testing whether a new study technique improves exam scores. The technique really works, but you only test a handful of students. You might not find convincing evidence of the improvement.'),
             div(class='explain-callout',strong('The key idea: '),'A study can miss a real effect. Statistical power is the probability that a study will produce a statistically significant result when a specific real effect exists.'),
             p(class='section-copy','PowerLab helps you understand how decisions such as the number of participants influence this probability.'),
             div(class='nav-actions',actionButton('intro_back1','← Back',class='btn-secondary'),actionButton('intro_next2','Learn the key concepts →',class='btn-start')))),
       tabPanel(title='Key concepts',value='concepts',
         div(class='stepper',span(class='step-pill','1 · Why power?'),span(class='step-pill active','2 · Key concepts'),span(class='step-pill','3 · Tests'),span(class='step-pill','4 · Try it'),span(class='step-pill','5 · Your level')),
         div(class='cardx',h2('The language of power analysis'),
             p(class='section-copy','Every study starts with a question. These six concepts help you understand how researchers turn that question into a statistical test. Explore the cards in any order; each has a more detailed explanation.'),
             div(class='learn-grid',
                 div(class='learn-card',h4('Sample size (n)'),p('How many independent observations or participants are included in your study.'),actionButton('concept_n','Learn more ⓘ',class='btn-secondary')),
                 div(class='learn-card',h4('Hypotheses (H₀ and H₁)'),p('H₀ describes a model without the tested effect; H₁ describes an effect or difference.'),actionButton('concept_hypothesis','Learn more ⓘ',class='btn-secondary')),
                 div(class='learn-card',h4('Significance and p-values'),p('A statistical test evaluates how unusual your data would be under H₀. A statistically significant result is not necessarily important in practice.'),actionButton('concept_pvalue','Learn more ⓘ',class='btn-secondary')),
                 div(class='learn-card',h4('Alpha (α)'),p('The chosen false-positive error rate for rejecting H₀ when it is true, often set at .05.'),actionButton('concept_alpha','Learn more ⓘ',class='btn-secondary')),
                 div(class='learn-card',h4('Effect size (d or f)'),p('A standardized way to describe how large the effect is assumed to be. Different statistical tests use different measures.'),actionButton('concept_d','Learn more ⓘ',class='btn-secondary')),
                 div(class='learn-card',h4('Power (1 − β)'),p('The chance that the chosen test is significant if the specific effect you assumed really exists.'),actionButton('concept_power','Learn more ⓘ',class='btn-secondary'))),
             div(class='nav-actions',actionButton('concepts_back','← Back',class='btn-secondary'),actionButton('concepts_next','Discover the tests →',class='btn-start')))),
       tabPanel(title='Tests explained',value='tests',
         div(class='stepper',span(class='step-pill','1 · Why power?'),span(class='step-pill','2 · Key concepts'),span(class='step-pill active','3 · Tests'),span(class='step-pill','4 · Try it'),span(class='step-pill','5 · Your level')),
         div(class='cardx',h2('Which statistical test fits your study?'),
             p(class='section-copy','The right test depends on what you want to compare. Choose a card to open a clear example, learn when to use the test, and understand which study designs are supported. Repeated-measures designs are not included.'),
             div(class='choice-grid',
                 actionButton('learn_one',HTML('<span class="choice-name">One-sample t-test</span><span class="choice-description">Compare one group mean with a fixed reference value.</span>'),class='choice-card'),
                 actionButton('learn_two',HTML('<span class="choice-name">Independent two-sample t-test</span><span class="choice-description">Compare two independent group means.</span>'),class='choice-card'),
                 actionButton('learn_anova',HTML('<span class="choice-name">One-way ANOVA</span><span class="choice-description">Compare means across several independent groups.</span>'),class='choice-card'),
                 actionButton('learn_anova2',HTML('<span class="choice-name">Two-way ANOVA</span><span class="choice-description">Investigate two factors and their interaction.</span>'),class='choice-card')),
             div(class='nav-actions',actionButton('tests_back','← Back',class='btn-secondary'),actionButton('tests_next','Try changing sample size →',class='btn-start')))),
       tabPanel(title='Interactive preview',value='try',
         div(class='stepper',span(class='step-pill','1 · Why power?'),span(class='step-pill','2 · Key concepts'),span(class='step-pill','3 · Tests'),span(class='step-pill active','4 · Try it'),span(class='step-pill','5 · Your level')),
         div(class='cardx',h2('Can a bigger study detect effects more reliably?'),
             p(class='section-copy','Imagine comparing exam results from students using two different study methods, and assume the new method truly raises scores by half a standard deviation (Cohen’s d = 0.5). The groups are compared with a two-sided independent t-test at α = .05. Move the slider to add more students to each group; the assumed effect and the test settings stay the same.'),
             div(class='intro-preview-slider',sliderInput('preview_n','Participants in each group',min=5,max=200,value=20,step=5,width='100%')),
             div(class='metric',textOutput('preview_power')),
             p(class='helper','This percentage is the chance that the test reports a significant result if the assumed difference really exists and the test assumptions hold.'),
             plotOutput('preview_curve',height='260px'),
             div(class='explain-callout',strong('What to notice: '),'For a fixed effect in the direction being tested, increasing the number of independent participants generally increases statistical power.'),
             div(class='nav-actions',actionButton('intro_back2','← Back',class='btn-secondary'),actionButton('intro_next3','Choose your learning level →',class='btn-start')))),
       tabPanel(title='Choose level',value='level',
         div(class='stepper',span(class='step-pill','1 · Why power?'),span(class='step-pill','2 · Key concepts'),span(class='step-pill','3 · Tests'),span(class='step-pill','4 · Try it'),span(class='step-pill active','5 · Your level')),
         div(class='cardx',h2('Choose your learning level'),
             p(class='section-copy','Choose Beginner for a guided introduction or Advanced to explore distributions, error probabilities and detailed statistics. Both modes use the same statistical calculations, and you can switch at any time.'),
             uiOutput('level_cards'),
             div(class='nav-actions',actionButton('intro_back3','← Back',class='btn-secondary'),actionButton('intro_finish','Enter PowerLab →',class='btn-start')))),
              # Explore page: one-screen dashboard. Parameters on the left, results and charts on the right,
              # so every change is visible immediately. Each slider has a linked number box (see slider_box()).
              tabPanel(title = 'Explore', value = 'explore',
                div(class='lab-page',
                       div(class='lab-hero',div(class='lab-hero-content',
                         div(div(class='lab-eyebrow','YOUR INTERACTIVE STATISTICS LAB'),
                             div(class='brand',HTML('▥ Power<span>Lab</span>'))),
                         div(class='lab-hero-actions',uiOutput('mode_switch'),
                             actionButton('home','← Introduction',class='btn-secondary')))),
                       div(style='display:none',selectInput('mode',NULL,choices=c('Beginner'='beginner','Advanced'='advanced'),selected='beginner')),
                       div(class='lab-main',
                         div(class='lab-left',
                           div(class='cardx lab-testcard',
                               card_title('① Choose your test','test'),
                               uiOutput('test_cards'),
                               div(style='display:none',radioButtons('test',NULL,c('One-sample t-test'='one','Independent two-sample t-test'='two', 'One-way ANOVA (F test)'='anova', 'Two-way ANOVA (F test)'='anova2'),selected='two')),
                               conditionalPanel(condition="input.test != 'anova' && input.test != 'anova2'",
                                                info_button('Alternative hypothesis','alternative'),
                                                radioButtons('alternative',NULL,c('Two-sided'='two','Greater'='greater','Less'='less'),selected='two',inline=TRUE))),
                           div(class='cardx lab-params',h3('② Set parameters'),
                               info_button('Sample size (n)','n'),
                               slider_box('n',3,200,30,1),
                               conditionalPanel(condition="input.test == 'anova'",
                                                info_button('Number of groups (k)','k'),
                                                slider_box('k',3,10,3,1)),
                               conditionalPanel(condition="input.test == 'anova2'",
                                                div(class='param-pair',
                                                    div(info_button('Levels of A','levels_a'),slider_box('levels_a',2,6,2,1)),
                                                    div(info_button('Levels of B','levels_b'),slider_box('levels_b',2,6,2,1))),
                                                info_button('Effect to test','effect'),
                                                radioButtons('effect',NULL,c('Main A'='A','Main B'='B','A × B'='AB'),selected='A',inline=TRUE)),
                               uiOutput('effect_label'),
                               slider_box('d',-2,2,0.5,.01),
                               info_button('Significance level (α)','alpha'),
                               slider_box('alpha',.001,.10,.05,.001),
                               info_button('Target power','target'),
                               slider_box('target',.7,.95,.8,.01),
                               uiOutput('direction_warning'))),
                         div(class='lab-right',
                           div(class='cardx lab-results',
                               div(class='card-title',h3('③ Your results'),
                                   actionButton('help_interpretation','What does this mean? ⓘ',class='btn-explain')),
                               div(class='metric-grid',
                                   div(info_button('Statistical power','power'),div(class='metric',textOutput('power',inline=TRUE))),
                                   div(info_button('Total participants','total'),div(class='metric',textOutput('total',inline=TRUE))),
                                   div(info_button('Required n','required'),div(class='metric',textOutput('required',inline=TRUE)))),
                               conditionalPanel(condition="input.mode != 'beginner'",
                                   div(class='lab-details',
                                       div(info_button('Degrees of freedom','df'),textOutput('df')),
                                       div(info_button('Noncentrality parameter','ncp'),textOutput('ncp')),
                                       div(info_button('Critical value(s)','critical'),textOutput('critical'))))),
                           div(class='lab-plots',
                             conditionalPanel(condition="input.mode != 'beginner'",class='dist-stack',
                               div(class='cardx plot-card',card_title('Explore the distributions','dist_note'),
                                   div(class='chart-key',actionButton('key_null',HTML('<span style="color:#58b8ff">●</span> H₀')),
                                       actionButton('key_alt',HTML('<span style="color:#ff729d">●</span> H₁')),
                                       actionButton('key_boundary','┊ Critical'),
                                       actionButton('key_region',HTML('<span style="color:#976cff">■</span> Power'))),
                                   plotOutput('distribution',height='260px',click=clickOpts(id='dist_click',clip=TRUE))),
                               div(class='cardx plot-card',card_title('Understand α, β and power','errors_note'),
                                   radioButtons('error_region',NULL,choices=c('Power (1 − β)'='power','Type II (β)'='beta','Type I (α)'='alpha','Correct (1 − α)'='correct'),selected='power',inline=TRUE),
                                   plotOutput('error_plot',height='260px'),div(class='error-caption',textOutput('error_caption')))),
                             div(class='curve-cell',
                               div(class='cardx plot-card',card_title('Power playground','curve_note'),
                                   div(class='chart-key',actionButton('key_curve',HTML('<span style="color:#a288ff">↗</span> Power curve')),
                                       actionButton('key_targetline',HTML('<span style="color:#ffb776">┄</span> Target')),
                                       actionButton('key_current',HTML('<span style="color:#a288ff">●</span> Your design'))),
                                   plotOutput('curve',height='340px',click=clickOpts(id='curve_click',clip=TRUE))))))))
              )
  )
)

server <- function(input, output, session) {
  # Test-specific teaching content. Designs here always involve independent observations.
  test_info <- list(
    one=list(title='One-sample t-test',detail='Compare one mean with a fixed value. Example: do students sleep more or less than 8 hours a night?'),
    two=list(title='Independent two-sample t-test',detail='Compare the means of two independent groups. Example: do students using a new study method score differently from those using the usual method?'),
    anova=list(title='One-way ANOVA',detail='Compare the means of three or more independent groups. Example: do exam scores differ between three teaching methods? A significant result does not tell you which groups differ.'),
    anova2=list(title='Two-way ANOVA',detail='Study two factors at once, and whether they interact. Example: do teaching method and class size affect exam scores, and does the effect of method depend on class size?'))
  level_choice <- reactiveVal('beginner')
  # Beginner and Advanced share all computations; only presentation differs.
  output$level_cards <- renderUI({
    div(class='choice-grid',lapply(list(
      list(id='beginner',title='🌱 Beginner',description='Learn the essentials with guided explanations and a clear power curve.'),
      list(id='advanced',title='🔬 Advanced',description='Explore sampling distributions, error probabilities, degrees of freedom and critical values.')
    ),function(z) actionButton(paste0('select_',z$id),
      HTML(paste0('<span class="choice-name">',z$title,'</span><span class="choice-description">',z$description,'</span>')),
      class=paste('choice-card',if(level_choice()==z$id)'selected' else ''))))
  })
  for (lev in c('beginner','advanced')) local({l<-lev;observeEvent(input[[paste0('select_',l)]],level_choice(l))})
  output$mode_switch <- renderUI({
    div(class='level-switch',role='group',`aria-label`='Learning level',
      actionButton('mode_beginner','Beginner',class=if(input$mode=='beginner')'active' else ''),
      actionButton('mode_advanced','Advanced',class=if(input$mode=='advanced')'active' else ''))
  })
  for (lev in c('beginner','advanced')) local({l<-lev;observeEvent(input[[paste0('mode_',l)]],{
    updateSelectInput(session,'mode',selected=l)
  })})
  output$test_cards <- renderUI({
    div(class='test-choice-list',lapply(names(test_info),function(id){
      z<-test_info[[id]]
      div(class='test-select-row',
          actionButton(paste0('choose_',id),HTML(paste0('<span class="choice-name">',z$title,'</span>')),
                       class=paste('choice-card test-card',if(input$test==id)'selected' else '')),
          actionButton(paste0('info_test_',id),'ⓘ',class='info test-info-button',title=paste('Learn about',z$title)))
    }))
  })
  for (tt in names(test_info)) local({t<-tt;observeEvent(input[[paste0('choose_',t)]],updateRadioButtons(session,'test',selected=t));observeEvent(input[[paste0('learn_',t)]],showModal(modalDialog(title=test_info[[t]]$title,p(test_info[[t]]$detail),easyClose=TRUE,footer=div(class='test-modal-actions',modalButton('Close'),actionButton(paste0('use_',t),'Choose this test',class='btn-primary')))));observeEvent(input[[paste0('use_',t)]],{updateRadioButtons(session,'test',selected=t);removeModal()})})
  for (tt in names(test_info)) local({t<-tt;observeEvent(input[[paste0('info_test_',t)]],showModal(modalDialog(title=test_info[[t]]$title,p(test_info[[t]]$detail),easyClose=TRUE,footer=modalButton('Close'))))})
  output$effect_label <- renderUI(info_button(if(is_f(input$test))'Effect size (Cohen’s f)' else 'Effect size (Cohen’s d)','d'))
  observeEvent(input$concepts_back,updateTabsetPanel(session,'page',selected='why'))
  observeEvent(input$concepts_next,updateTabsetPanel(session,'page',selected='tests'))
  observeEvent(input$tests_back,updateTabsetPanel(session,'page',selected='concepts'))
  observeEvent(input$tests_next,updateTabsetPanel(session,'page',selected='try'))
  for (k in c('n','hypothesis','pvalue','alpha','d','power')) local({key<-k;observeEvent(input[[paste0('concept_',key)]],{z<-help_text[[key]];showModal(modalDialog(title=z[[1]],p(z[[2]]),easyClose=TRUE,footer=modalButton('Got it')))})})

  # Switch the effect-size slider between Cohen's d (t-tests) and Cohen's f (ANOVA)
  was_f <- reactiveVal(FALSE)
  observeEvent(input$test, {
    is_a <- is_f(input$test)
    if (is_a != was_f()) {
      if (is_a) updateSliderInput(session,'d',min=0,max=1,value=0.25,step=.01)
      else updateSliderInput(session,'d',min=-2,max=2,value=0.5,step=.01)
      # Update the box's value together with its range; a range-only update would make the box
      # re-send its old value, which would then override the slider.
      updateNumericInput(session,'d_num',value=if(is_a)0.25 else 0.5,min=if(is_a)0 else -2,max=if(is_a)1 else 2)
      was_f(is_a)
    }
  })

  # Keep each slider and its number box in sync. Typing is debounced so a value like 120 is not
  # applied digit by digit; typed values are snapped to the slider's range and step.
  for (id in names(slider_limits)) local({
    s <- id; b <- paste0(id,'_num')
    typed <- debounce(reactive(input[[b]]), 600)
    observeEvent(typed(), {
      v <- typed(); req(is.numeric(v), length(v)==1, !is.na(v))
      lim <- slider_limits[[s]]
      if (s=='d' && is_f(input$test)) lim[1:2] <- c(0,1)
      v <- min(max(v,lim[1]),lim[2])
      v <- round(round((v-lim[1])/lim[3])*lim[3]+lim[1],6)
      if (!isTRUE(all.equal(v,input[[s]]))) updateSliderInput(session,s,value=v)
      else if (!isTRUE(all.equal(v,typed()))) updateNumericInput(session,b,value=v)
    }, ignoreInit=TRUE)
    observeEvent(input[[s]], {
      cur <- isolate(input[[b]])
      if (is.null(cur) || is.na(cur) || !isTRUE(all.equal(cur,input[[s]]))) updateNumericInput(session,b,value=round(input[[s]],6))
    })
  })
  
  observeEvent(input$intro_next1,updateTabsetPanel(session,'page',selected='why'))
  observeEvent(input$intro_next2,updateTabsetPanel(session,'page',selected='concepts'))
  observeEvent(input$intro_next3,updateTabsetPanel(session,'page',selected='level'))
  observeEvent(input$intro_back1,updateTabsetPanel(session,'page',selected='intro'))
  observeEvent(input$intro_back2,updateTabsetPanel(session,'page',selected='tests'))
  observeEvent(input$intro_back3,updateTabsetPanel(session,'page',selected='try'))
  observeEvent(input$intro_skip,{updateSelectInput(session,'mode',selected='beginner');updateTabsetPanel(session,'page',selected='explore')})
  observeEvent(input$intro_finish,{updateSelectInput(session,'mode',selected=level_choice());updateTabsetPanel(session,'page',selected='explore')})
  output$preview_power <- renderText(sprintf('%.1f%% power',100*calc_power(input$preview_n,.5,.05,'two','two')$power))
  output$preview_curve <- renderPlot({dd <- data.frame(n=5:200,power=vapply(5:200,function(n)calc_power(n,.5,.05,'two','two')$power,numeric(1)));ggplot(dd,aes(n,power))+geom_line(color='#aa8dff',linewidth=1.3)+geom_point(data=subset(dd,n==input$preview_n),color='#ff729d',size=4)+labs(x='Participants per group',y='Power')+theme_minimal(base_size=15)+theme(text=element_text(size=15,color='#dce5f9'),axis.text=element_text(size=13,color='#dce5f9'),plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))},bg='#121b32')
  observeEvent(input$home,updateTabsetPanel(session,'page',selected='intro'))
  show_explanation <- function(key, extra = NULL) {
    item <- help_text[[key]]
    if (is.null(item)) return(invisible(NULL))
    showModal(modalDialog(title=item[[1]],p(item[[2]]),
                          if (!is.null(extra)) p(tags$strong(extra)),
                          easyClose=TRUE,footer=modalButton('Got it')))
  }
  for (key in c('test','alternative','n','k','levels_a','levels_b','effect','d','alpha','target','power','total','required','df','ncp','critical','dist_note','errors_note','curve_note')) {
    local({ k <- key; observeEvent(input[[paste0('help_',k)]],show_explanation(k)) })
  }
  for (key in c('null','alt','boundary','region','curve','targetline','current')) {
    local({ k <- key; observeEvent(input[[paste0('key_',k)]],show_explanation(k)) })
  }
  # The plain-language interpretation now opens in a pop-up; it updates live if parameters change.
  observeEvent(input$help_interpretation,
    showModal(modalDialog(title='What do your results mean?',textOutput('interpretation'),
                          easyClose=TRUE,footer=modalButton('Got it'))))
  
  # Extra design information passed to calc_power(): k for one-way, list(a, b, effect) for two-way
  design <- reactive(if (input$test == 'anova2') list(a = input$levels_a, b = input$levels_b, effect = input$effect) else input$k)
  # Wording for the sampling unit
  unit <- reactive(if (input$test == 'anova2') 'cell' else 'group')
  
  result <- reactive(calc_power(input$n,input$d,input$alpha,input$test,input$alternative,design()))
  output$power <- renderText(sprintf('%.1f%%',100*result()$power))
  output$total <- renderText(result()$total)
  output$df <- renderText(if(is_f(input$test)) sprintf('df₁ = %d, df₂ = %d',result()$df1,result()$df) else result()$df)
  output$ncp <- renderText(sprintf(if(is_f(input$test)) 'λ = %.3f' else 'δ = %.3f',result()$delta))
  output$critical <- renderText(paste(sprintf('%+.3f',result()$critical),collapse='  /  '))
  
  # Search up to n = 5000 to avoid silently treating 200 as a universal cap.
  recommendation <- reactive({
    if (!is_f(input$test) && ((input$alternative=='greater' && input$d<0) || (input$alternative=='less' && input$d>0))) return(NA_integer_)
    # pt() warns about reduced precision for very large noncentrality values; results remain accurate here.
    vals <- suppressWarnings(vapply(3:5000,function(n) calc_power(n,input$d,input$alpha,input$test,input$alternative,design())$power,numeric(1)))
    hit <- which(vals >= input$target)
    if (length(hit)==0) NA_integer_ else hit[1]+2L
  })
  direction_conflict <- reactive(!is_f(input$test) && ((input$alternative=='greater' && input$d<0) || (input$alternative=='less' && input$d>0)))
  output$direction_warning <- renderUI(if(direction_conflict())div(class='warning-note','The effect size sign and one-sided alternative conflict. The actual power is shown, but no sample-size recommendation is meaningful until you correct the direction.') else NULL)
  output$required <- renderText({
    if (is.na(recommendation()) && direction_conflict()) 'Direction mismatch' else if (is.na(recommendation())) '> 5,000' else if (input$test!='one') paste0(recommendation(),' / ',unit()) else as.character(recommendation())
  })
  output$interpretation <- renderText({
    r <- result(); req <- recommendation()
    eff <- if(is_f(input$test)) 'f' else 'd'
    mult <- switch(input$test,one=1,two=2,anova=input$k,anova2=input$levels_a*input$levels_b)
    recommendation_msg <- if(direction_conflict()) 'Warning: the assumed effect points opposite to the selected one-sided hypothesis. Increasing n will not solve this directional mismatch; check the hypothesis and effect sign.' else if(is.na(req)) 'The target is not reached for n up to 5,000 per sampling unit in the tested direction.' else
      if(input$test!='one') sprintf('To reach the target, plan for at least %d per %s (%d participants total).',req,unit(),mult*req) else
        sprintf('To reach the target, plan for at least %d observations.',req)
    paste0(sprintf('If the true effect is %s = %.2f, this design has %.1f%% power for the selected test. ',eff,input$d,r$power*100),
           if(r$power >= input$target) 'Your selected target is met. ' else 'Your selected target is not yet met. ',
           recommendation_msg)
  })
  
  dist_data <- reactive({
    r <- result(); tst <- input$test
    xx <- if (is_f(tst)) seq(0.001, max(qf(.999,r$df1,r$df), r$critical*1.2, qf(.999,r$df1,r$df,ncp=r$delta)), length.out=2401)
    else seq(min(-6,r$delta-6,min(r$critical)-1),max(6,r$delta+6,max(r$critical)+1),length.out=2401)
    h1 <- dens1(xx,r,tst)
    data.frame(x=xx,h0=dens0(xx,r,tst),h1=h1,
               shaded=ifelse(reject_region(xx,r,input$alternative,tst),h1,0))
  })
  output$distribution <- renderPlot({
    dat <- dist_data(); r <- result()
    ggplot(dat,aes(x=x))+
      geom_area(aes(y=shaded),fill='#976cff',alpha=.48)+
      geom_line(aes(y=h0),color='#58b8ff',linewidth=1.25)+
      geom_line(aes(y=h1),color='#ff729d',linewidth=1.25)+
      geom_vline(xintercept=r$critical,color='#c7d1ea',linetype='dashed',linewidth=.7)+
      labs(x=stat_label(input$test),y='Density')+
      coord_cartesian(ylim=c(0,y_top(dat,r,input$test)))+
      theme_minimal(base_size=15)+
      theme(text=element_text(size=15,color='#dce5f9'),axis.text=element_text(size=13,color='#aebfdf'),
            panel.grid=element_line(color='#2b3853'),panel.grid.minor=element_blank(),
            plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))
  },bg='#121b32')
  observeEvent(input$dist_click,{
    click <- input$dist_click; dat <- dist_data(); r <- result()
    if(is.null(click$x)||is.null(click$y)) return()
    yrange <- y_top(dat,r,input$test)
    xr <- range(dat$x); x <- click$x; y <- click$y
    # Prioritize close lines first, then critical boundaries, then shaded region.
    h0 <- dens0(x,r,input$test); h1 <- dens1(x,r,input$test)
    d0 <- abs(y-h0)/yrange; d1 <- abs(y-h1)/yrange
    crit_distance <- min(abs(x-r$critical))/diff(xr)
    if(min(d0,d1)<.055) show_explanation(if(d0<=d1)'null' else 'alt')
    else if(crit_distance<.025) show_explanation('boundary')
    else if(y>0 && y<h1 && reject_region(x,r,input$alternative,input$test)) show_explanation('region')
  })
  
  # Show n = 3..200 by default; extend the axis when the required n lies beyond 200,
  # so the curve always visibly crosses the target line when the target is reachable.
  curve_max <- reactive({
    req_n <- recommendation()
    if (is.na(req_n) || req_n <= 180) 200L else as.integer(min(5000, ceiling(req_n * 1.15)))
  })
  curve_data <- reactive({
    nn <- 3:curve_max()
    data.frame(n=nn,power=suppressWarnings(vapply(nn,function(n)calc_power(n,input$d,input$alpha,input$test,input$alternative,design())$power,numeric(1))))
  })
  output$curve <- renderPlot({
    dat <- curve_data()
    ggplot(dat,aes(n,power))+
      geom_hline(yintercept=input$target,color='#ffb776',linetype='dashed',linewidth=.75)+
      geom_line(color='#a288ff',linewidth=1.35)+
      geom_point(data=data.frame(n=input$n,power=result()$power),color='#ffffff',fill='#a288ff',shape=21,size=4,stroke=1.2)+
      scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(x*100),'%'))+
      labs(x=if(input$test=='one')'Total sample size' else paste('Participants per',unit()),y='Power')+
      theme_minimal(base_size=15)+
      theme(text=element_text(size=15,color='#dce5f9'),axis.text=element_text(size=13,color='#aebfdf'),
            panel.grid=element_line(color='#2b3853'),panel.grid.minor=element_blank(),
            plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))
  },bg='#121b32')
  output$error_plot <- renderPlot({
    dat <- dist_data(); r <- result(); a <- input$error_region
    reject <- reject_region(dat$x,r,input$alternative,input$test)
    dat$region <- if(a=='alpha') ifelse(reject,dat$h0,0) else if(a=='correct') ifelse(!reject,dat$h0,0) else if(a=='beta') ifelse(!reject,dat$h1,0) else ifelse(reject,dat$h1,0)
    ggplot(dat,aes(x))+geom_area(aes(y=region),fill='#9c79ff',alpha=.55)+geom_line(aes(y=h0),color='#58b8ff',linewidth=1)+geom_line(aes(y=h1),color='#ff729d',linewidth=1)+geom_vline(xintercept=r$critical,linetype='dashed',color='#e0e5ff')+labs(x=stat_label(input$test),y='Density')+coord_cartesian(ylim=c(0,y_top(dat,r,input$test)))+theme_minimal(base_size=15)+theme(text=element_text(size=15,color='#dce5f9'),axis.text=element_text(size=13,color='#dce5f9'),panel.grid=element_line(color='#2b3853'),plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))
  },bg='#121b32')
  output$error_caption <- renderText({r <- result();switch(input$error_region,alpha=sprintf('Type I error: %.1f%% of null-distribution area is in the rejection region.',100*input$alpha),correct=sprintf('No Type I error under H₀: %.1f%%.',100*(1-input$alpha)),beta=sprintf('Type II error at the specified effect: %.1f%%.',100*(1-r$power)),power=sprintf('Power at the specified effect: %.1f%%.',100*r$power))})
  observeEvent(input$curve_click,{
    click <- input$curve_click
    if(is.null(click$x)||is.null(click$y)) return()
    nmax <- curve_max()
    if(abs(click$x-input$n)<max(5,.025*nmax) && abs(click$y-result()$power)<.08) show_explanation('current')
    else if(abs(click$y-input$target)<.035) show_explanation('targetline')
    else {
      closest_n <- max(3,min(nmax,round(click$x)))
      expected <- calc_power(closest_n,input$d,input$alpha,input$test,input$alternative,design())$power
      if(abs(click$y-expected)<.09) show_explanation('curve',sprintf('At n = %d, predicted power = %.1f%%.',closest_n,100*expected))
    }
  })
}
shinyApp(ui,server)
