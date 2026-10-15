library(shiny)
library(ggplot2)

# PowerLab: the current version covers one-sample and equal-size,
# equal-variance independent two-sample t tests.
calc_power <- function(n, d, alpha, test, alternative) {
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

help_text <- list(
  test = c('Which t-test?', 'A one-sample t-test compares a sample mean with a fixed reference value. An independent two-sample t-test compares the means of two independent groups. This app uses the pooled-variance independent t-test with equal group sizes and assumes equal population variances.'),
  alternative = c('Alternative hypothesis', 'Two-sided tests can reject for sufficiently large differences in either direction. Greater tests reject for positive standardized differences; less tests reject for negative ones. Choose the direction before examining the data.'),
  n = c('Sample size', 'For one-sample tests, n is the number of observations. For independent two-sample tests, n is the number in EACH group (total N = 2n). Larger samples generally increase power for a fixed nonzero effect.'),
  d = c("Cohen’s d", 'The true standardized difference assumed for planning: (population mean − reference)/population SD for one sample; (population mean 1 − population mean 2)/common population SD for two groups. Negative d indicates the opposite direction. Values such as 0.2/0.5/0.8 are conventions, not universal scientific thresholds.'),
  alpha = c('Significance level (α)', 'Alpha is the long-run probability of rejecting the null hypothesis when it is true (a Type I error), under the test assumptions. It is not the probability that H₀ is true. Smaller α makes rejection harder.'),
  target = c('Target power', 'Your desired power at the specified true effect. An 80% target is common, not a universal rule. The recommended N is the smallest integer sample size searched that reaches the target, under the assumptions you entered.'),
  power = c('Statistical power (1 − β)', 'The probability of rejecting H₀ if the SPECIFIED population effect is true. For example, 80% power means approximately 80 of 100 independent repetitions would reject H₀ if the planning assumptions are correct. It is not the probability that your hypothesis is correct.'),
  total = c('Total sample size', 'For the one-sample t-test, total N = n. For an independent two-sample t-test in this app, the groups are equally sized, so total N = 2n.'),
  required = c('Required sample size', 'This is the smallest integer n that reaches the chosen target power in the search range. For a two-sample test it is n PER GROUP. If a target cannot be reached within the search limit, that is reported explicitly.'),
  df = c('Degrees of freedom', 'One-sample t-test: df = n − 1. Equal-size pooled two-sample t-test: df = 2n − 2. Degrees of freedom affect the shape and critical values of the t-distribution.'),
  ncp = c('Noncentrality parameter (δ)', 'Under the specified alternative, the t statistic follows a noncentral t-distribution with δ = d√n (one sample) or δ = d√(n/2) (two equal independent groups). This is not an effect size estimate.'),
  critical = c('Critical t-value', 'The t-statistic threshold beyond which a result would be statistically significant. A two-sided test has two thresholds, −tcrit and +tcrit; a one-sided test has one. The threshold depends on α, df and test direction.'),
  null = c('Blue curve: H₀', 'The central t-distribution of the t-statistic under the null hypothesis of zero population difference. Its area across all t-values is 1; it represents repeated-sample variability of the statistic, not raw participant scores.'),
  alt = c('Pink curve: H₁', 'The noncentral t-distribution of the t-statistic under the specific true effect (Cohen’s d) currently entered. The distribution shifts according to the noncentrality parameter and determines statistical power.'),
  region = c('Purple shading: power', 'The shaded area is the probability that a t-statistic falls inside the rejection region when the specified alternative is true. That probability is statistical power (1 − β). For a two-sided test, BOTH rejection tails contribute.'),
  boundary = c('Dashed line: rejection threshold', 'The vertical dashed line marks a critical t-value. Values beyond the relevant boundary (or beyond either boundary for a two-sided test) lead to rejection of H₀. Boundaries are based on the null distribution.'),
  curve = c('Power vs sample size', 'The curve recalculates power for different n, holding Cohen’s d, α, test type and alternative direction constant. The highlighted point is your current design. In a two-group test the horizontal axis is n PER GROUP.'),
  targetline = c('Dashed horizontal line: target', 'The dashed horizontal line represents your chosen target power. Where the power curve first reaches it gives the smallest sample size in the searched range that meets the target.'),
  type1 = c('Type I error (α)', 'The probability of rejecting H₀ when it is true. The area of the null distribution in the rejection region equals α.'),
  type2 = c('Type II error (β)', 'The probability of failing to reject H₀ when the specified alternative is true. It equals 1 minus power.'),
  sensitivity = c('Sensitivity analysis', 'With sample size, alpha and desired power fixed, sensitivity analysis finds the smallest effect magnitude detectable with that target power, for a prespecified direction.'),
  current = c('Highlighted point: your design', 'The marker shows your current sample size and its calculated power. Adjust the n slider to watch the marker move along the power curve.')
)

info_button <- function(text, key) {
  div(class = 'labelrow', tags$label(text),
      actionButton(paste0('help_', key), 'ⓘ', class = 'info', title = paste('Learn about', text)))
}

ui <- fluidPage(
  tags$head(
    tags$meta(name = 'viewport', content = 'width=device-width, initial-scale=1'),
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
      @media(prefers-reduced-motion:reduce){*,*:before,*:after{animation:none!important;transition:none!important}}
    "))
  ),
  div(class = 'topbar',
      div(div(class = 'brand', HTML('▥ Power<span>Lab</span>')),div(class = 'subbrand','Explore · Understand · Plan')),
      div(class = 'steps',span(class='step','Learn'),span(class='step','Explore'),span(class='step','Apply'))),
  tabsetPanel(id = 'page', type = 'hidden',
    tabPanel(title='Introduction',value='intro',
      div(class='hero',div(class='eyebrow','Welcome to PowerLab'),
          h1('Statistics you can ',span('see.')),
          p('Learn how study design choices influence your probability of detecting an effect. Experiment with the statistical distributions, rather than memorizing formulas.'),
          div(class='intro-graphic',lapply(c(24,42,67,90,67,42,24),function(h) tags$span(style=paste0('height:',h,'px')))),
          actionButton('intro_next1','Discover why power matters →',class='btn-start'),
          actionButton('intro_skip','Skip introduction',class='btn-secondary'))),
    tabPanel(title='Why power?',value='why',
      div(class='cardx',div(class='eyebrow','02 / 04 · Motivation'),h2('Why does statistical power matter?'),
          p('Imagine testing whether a stress-reduction intervention improves wellbeing. The intervention really has a positive effect, but your study uses very few participants. Your statistical test may fail to detect that effect.'),
          p('A non-significant result does not establish that there is no effect. Power analysis helps you plan the sample size needed to detect a specific effect with a chosen probability.'),
          div(class='mini',strong('Think about it'),p('What would happen if you recruit twice as many independent participants, holding the effect size and alpha constant?'),
              actionButton('intro_reveal','Reveal the answer',class='btn-secondary'),uiOutput('intro_answer')),
          actionButton('intro_back1','← Back',class='btn-secondary'),actionButton('intro_next2','Try it yourself →',class='btn-start'))),
    tabPanel(title='Try it',value='try',
      div(class='cardx',div(class='eyebrow','03 / 04 · Interactive preview'),h2('What changes statistical power?'),
          p('Move the sample-size slider. The curve represents statistical power for a two-sided independent-samples t-test with a true standardized effect d = 0.5, α = .05, and equally sized groups.'),
          sliderInput('preview_n','Participants per group',min=5,max=200,value=20,step=5),
          div(class='metric',textOutput('preview_power')),p(class='helper','This is the probability of a significant result if the assumed alternative and model are true.'),
          plotOutput('preview_curve',height='220px'),
          actionButton('intro_back2','← Back',class='btn-secondary'),actionButton('intro_next3','Choose your learning level →',class='btn-start'))),
    tabPanel(title='Choose level',value='level',
      div(class='cardx',div(class='eyebrow','04 / 04 · Your experience'),h2('Choose your level'),
          p('You can change the level at any time. The statistics stay the same; the level changes how much explanation and how many tools are shown.'),
          fluidRow(
            column(4,div(class='mini',h3('Beginner'),p('Guided terminology, essential inputs, and approachable explanations.'))),
            column(4,div(class='mini',h3('Intermediate'),p('Explore tails, error regions, comparisons, and underlying statistics.'))),
            column(4,div(class='mini',h3('Expert'),p('Exact planning parameters, sensitivity analysis, and reproducible exports.')))),
          radioButtons('intro_mode','Your starting level',c('Beginner'='beginner','Intermediate'='intermediate','Expert'='expert'),selected='beginner',inline=TRUE),
          actionButton('intro_back3','← Back',class='btn-secondary'),actionButton('intro_finish','Start PowerLab →',class='btn-start'))),
    tabPanel(title = 'Explore', value = 'explore',
      fluidRow(column(12,div(class='cardx',fluidRow(column(5,div(class='eyebrow','LEARN · EXPLORE · APPLY'),h3('Your interactive statistics lab')),column(4,selectInput('mode','Learning level',choices=c('Beginner'='beginner','Intermediate'='intermediate','Expert'='expert'),selected='beginner')),column(3,actionButton('home','← Introduction',class='btn-secondary'))),textOutput('mode_description')))),
      fluidRow(
        column(4,
          div(class='cardx',
              h3(style='margin-top:20px','① Choose your test'),
              info_button('Statistical test','test'),
              radioButtons('test',NULL,c('One-sample t-test'='one','Independent two-sample t-test'='two'),selected='two'),
              info_button('Alternative hypothesis','alternative'),
              radioButtons('alternative',NULL,c('Two-sided'='two','Greater (positive)'='greater','Less (negative)'='less'),selected='two')),
          div(class='cardx',h3('② Set parameters'),
              info_button('Sample size (n)','n'),
              sliderInput('n',NULL,3,200,30,step=1),
              info_button('Effect size (Cohen’s d)','d'),
              sliderInput('d',NULL,-2,2,0.5,step=.05),
              info_button('Significance level (α)','alpha'),
              sliderInput('alpha',NULL,.001,.10,.05,step=.001),
              info_button('Target power','target'),
              sliderInput('target',NULL,.7,.95,.8,step=.05))
        ),
        column(8,
          div(class='cardx',h3('③ Explore your results'),
              fluidRow(
                column(4,info_button('Statistical power','power'),div(class='metric',textOutput('power',inline=TRUE))),
                column(4,info_button('Total participants','total'),div(class='metric',textOutput('total',inline=TRUE))),
                column(4,info_button('Required n','required'),div(class='metric',textOutput('required',inline=TRUE))))),
          div(class='cardx',h3('Sampling distributions'),
              p(class='chart-note','Click a curve, dashed boundary, or shaded region to learn what it represents.'),
              div(class='chart-key',actionButton('key_null','● H₀ curve'),actionButton('key_alt','● H₁ curve'),actionButton('key_boundary','┊ Critical value'),actionButton('key_region','■ Power area')),
              plotOutput('distribution',height='295px',click=clickOpts(id='dist_click',clip=TRUE))),
          div(class='cardx',h3('Power changes with sample size'),
              p(class='chart-note','Click the curve or highlighted dot; the dashed horizontal line marks your target.'),
              div(class='chart-key',actionButton('key_curve','↗ Power curve'),actionButton('key_targetline','┄ Target power'),actionButton('key_current','● Your current design')),
              plotOutput('curve',height='250px',click=clickOpts(id='curve_click',clip=TRUE))),
          div(class='cardx',h3('Interpretation'),textOutput('interpretation'),hr(),
              conditionalPanel(condition="input.mode != 'beginner'",info_button('Degrees of freedom','df'),textOutput('df'),
              info_button('Noncentrality parameter','ncp'),textOutput('ncp'),
              info_button('Critical t-value(s)','critical'),textOutput('critical')))
        )
      ),
      conditionalPanel(condition="input.mode != 'beginner'",
        div(class='cardx',h3('Understand α, β and power'),p(class='helper','Select the highlighted region to understand what kind of decision it represents.'),
        radioButtons('error_region','Highlight',choices=c('Power (1 − β)'='power','Type II error (β)'='beta','Type I error (α)'='alpha','Correct non-rejection (1 − α)'='correct'),selected='power',inline=TRUE),
        plotOutput('error_plot',height='250px'),textOutput('error_caption'))),
      conditionalPanel(condition="input.mode != 'beginner'",div(class='cardx',h3('What happens if…?'),p(class='helper','Compare your current design with a changed sample size. Effect size, alpha and the test remain fixed.'),
        sliderInput('compare_n','Comparison: sample size',min=3,max=500,value=60,step=1),
        tableOutput('comparison'),actionButton('challenge','Try a learning challenge',class='btn-secondary'),textOutput('challenge_text'))),
      conditionalPanel(condition="input.mode == 'expert'",div(class='cardx',h3('Apply to your research'),
        p('Sensitivity analysis: the smallest positive effect detectable at your target power with the current n and α, for a two-sided or greater-direction test.'),
        textOutput('sensitivity_result'),downloadButton('download_report','Download analysis report'),downloadButton('download_csv','Download power curve CSV'),
        p(class='helper','The exported report records all selected assumptions.'))),
      p(class='footnote','Planning assumptions: independent observations; normally distributed observations within each group for the exact t calculations; two independent groups are equal size and have equal population variance. The effect size is assumed known for planning. Designed as an educational tool, not a substitute for checking study assumptions.')
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$intro_next1,updateTabsetPanel(session,'page',selected='why'))
  observeEvent(input$intro_next2,updateTabsetPanel(session,'page',selected='try'))
  observeEvent(input$intro_next3,updateTabsetPanel(session,'page',selected='level'))
  observeEvent(input$intro_back1,updateTabsetPanel(session,'page',selected='intro'))
  observeEvent(input$intro_back2,updateTabsetPanel(session,'page',selected='why'))
  observeEvent(input$intro_back3,updateTabsetPanel(session,'page',selected='try'))
  observeEvent(input$intro_skip,{updateSelectInput(session,'mode',selected='beginner');updateTabsetPanel(session,'page',selected='explore')})
  observeEvent(input$intro_finish,{updateSelectInput(session,'mode',selected=input$intro_mode);updateTabsetPanel(session,'page',selected='explore')})
  observeEvent(input$intro_reveal,output$intro_answer <- renderUI(p('With a fixed nonzero true effect and the same test assumptions, a larger sample generally means higher power.')))
  output$preview_power <- renderText(sprintf('%.1f%% power',100*calc_power(input$preview_n,.5,.05,'two','two')$power))
  output$preview_curve <- renderPlot({dd <- data.frame(n=5:200,power=vapply(5:200,function(n)calc_power(n,.5,.05,'two','two')$power,numeric(1)));ggplot(dd,aes(n,power))+geom_line(color='#aa8dff',linewidth=1.3)+geom_point(data=subset(dd,n==input$preview_n),color='#ff729d',size=4)+labs(x='Participants per group',y='Power')+theme_minimal()+theme(text=element_text(color='#dce5f9'),axis.text=element_text(color='#dce5f9'),plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))},bg='#121b32')
  output$mode_description <- renderText(switch(input$mode,beginner='Beginner · Focus on the key concepts and interactive plots.',intermediate='Intermediate · Explore statistical errors and compare sample sizes.',expert='Expert · Study error regions, compare designs and export reproducible analyses.'))
  observeEvent(input$home,updateTabsetPanel(session,'page',selected='intro'))
  show_explanation <- function(key, extra = NULL) {
    item <- help_text[[key]]
    if (is.null(item)) return(invisible(NULL))
    showModal(modalDialog(title=item[[1]],p(item[[2]]),
      if (!is.null(extra)) p(tags$strong(extra)),
      easyClose=TRUE,footer=modalButton('Got it')))
  }
  for (key in c('test','alternative','n','d','alpha','target','power','total','required','df','ncp','critical')) {
    local({ k <- key; observeEvent(input[[paste0('help_',k)]],show_explanation(k)) })
  }
  for (key in c('null','alt','boundary','region','curve','targetline','current')) {
    local({ k <- key; observeEvent(input[[paste0('key_',k)]],show_explanation(k)) })
  }
  result <- reactive(calc_power(input$n,input$d,input$alpha,input$test,input$alternative))
  output$power <- renderText(sprintf('%.1f%%',100*result()$power))
  output$total <- renderText(result()$total)
  output$df <- renderText(result()$df)
  output$ncp <- renderText(sprintf('δ = %.3f',result()$delta))
  output$critical <- renderText(paste(sprintf('%+.3f',result()$critical),collapse='  /  '))

  # Search up to n = 5000 to avoid silently treating 200 as a universal cap.
  recommendation <- reactive({
    vals <- vapply(3:5000,function(n) calc_power(n,input$d,input$alpha,input$test,input$alternative)$power,numeric(1))
    hit <- which(vals >= input$target)
    if (length(hit)==0) NA_integer_ else hit[1]+2L
  })
  output$required <- renderText({
    if (is.na(recommendation())) '> 5,000' else if (input$test=='two') paste0(recommendation(),' / group') else as.character(recommendation())
  })
  output$interpretation <- renderText({
    r <- result(); req <- recommendation()
    recommendation_msg <- if(is.na(req)) 'The target is not reached for n up to 5,000 in the tested direction.' else
      if(input$test=='two') sprintf('To reach the target, plan for at least %d per group (%d participants total).',req,2*req) else
        sprintf('To reach the target, plan for at least %d observations.',req)
    paste0(sprintf('If the true effect is d = %.2f, this design has %.1f%% power for the selected test. ',input$d,r$power*100),
           if(r$power >= input$target) 'Your selected target is met. ' else 'Your selected target is not yet met. ',
           recommendation_msg)
  })

  dist_data <- reactive({
    r <- result()
    xx <- seq(min(-6,r$delta-6,min(r$critical)-1),max(6,r$delta+6,max(r$critical)+1),length.out=2401)
    reject <- if(input$alternative=='two') abs(xx)>=max(abs(r$critical)) else if(input$alternative=='greater') xx>=r$critical else xx<=r$critical
    data.frame(x=xx,h0=dt(xx,df=r$df),h1=dt(xx,df=r$df,ncp=r$delta),shaded=ifelse(reject,dt(xx,df=r$df,ncp=r$delta),0))
  })
  output$distribution <- renderPlot({
    dat <- dist_data(); r <- result()
    ggplot(dat,aes(x=x))+
      geom_area(aes(y=shaded),fill='#976cff',alpha=.48)+
      geom_line(aes(y=h0),color='#58b8ff',linewidth=1.25)+
      geom_line(aes(y=h1),color='#ff729d',linewidth=1.25)+
      geom_vline(xintercept=r$critical,color='#c7d1ea',linetype='dashed',linewidth=.7)+
      labs(x='t-statistic',y='Density')+
      coord_cartesian(ylim=c(0,max(dat$h0,dat$h1)*1.12))+
      theme_minimal(base_size=12)+
      theme(text=element_text(color='#dce5f9'),axis.text=element_text(color='#aebfdf'),
            panel.grid=element_line(color='#2b3853'),panel.grid.minor=element_blank(),
            plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))
  },bg='#121b32')
  observeEvent(input$dist_click,{
    click <- input$dist_click; dat <- dist_data();r <- result()
    if(is.null(click$x)||is.null(click$y)) return()
    yrange <- max(dat$h0,dat$h1)*1.12
    xr <- range(dat$x); x <- click$x; y <- click$y
    # Prioritize close lines first, then critical boundaries, then shaded region.
    h0 <- dt(x,df=r$df); h1 <- dt(x,df=r$df,ncp=r$delta)
    d0 <- abs(y-h0)/yrange; d1 <- abs(y-h1)/yrange
    crit_distance <- min(abs(x-r$critical))/diff(xr)
    if(min(d0,d1)<.055) show_explanation(if(d0<=d1)'null' else 'alt')
    else if(crit_distance<.025) show_explanation('boundary')
    else if(y>0 && y<h1 && (if(input$alternative=='two') abs(x)>=abs(r$critical[2]) else if(input$alternative=='greater') x>=r$critical else x<=r$critical)) show_explanation('region')
  })

  curve_data <- reactive({
    nn <- 3:200
    data.frame(n=nn,power=vapply(nn,function(n)calc_power(n,input$d,input$alpha,input$test,input$alternative)$power,numeric(1)))
  })
  output$curve <- renderPlot({
    dat <- curve_data()
    ggplot(dat,aes(n,power))+
      geom_hline(yintercept=input$target,color='#ffb776',linetype='dashed',linewidth=.75)+
      geom_line(color='#a288ff',linewidth=1.35)+
      geom_point(data=data.frame(n=input$n,power=result()$power),color='#fff',fill='#a288ff',shape=21,size=4,stroke=1.2)+
      scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(x*100),'%'))+
      labs(x=if(input$test=='one')'Total sample size' else 'Sample size per group',y='Power')+
      theme_minimal(base_size=12)+
      theme(text=element_text(color='#dce5f9'),axis.text=element_text(color='#aebfdf'),
            panel.grid=element_line(color='#2b3853'),panel.grid.minor=element_blank(),
            plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))
  },bg='#121b32')
  output$error_plot <- renderPlot({
    dat <- dist_data(); r <- result(); a <- input$error_region
    reject <- if(input$alternative=='two') abs(dat$x)>=max(abs(r$critical)) else if(input$alternative=='greater') dat$x>=r$critical else dat$x<=r$critical
    dat$region <- if(a=='alpha') ifelse(reject,dat$h0,0) else if(a=='correct') ifelse(!reject,dat$h0,0) else if(a=='beta') ifelse(!reject,dat$h1,0) else ifelse(reject,dat$h1,0)
    ggplot(dat,aes(x))+geom_area(aes(y=region),fill='#9c79ff',alpha=.55)+geom_line(aes(y=h0),color='#58b8ff',linewidth=1)+geom_line(aes(y=h1),color='#ff729d',linewidth=1)+geom_vline(xintercept=r$critical,linetype='dashed',color='#e0e5ff')+labs(x='t-statistic',y='Density')+theme_minimal()+theme(text=element_text(color='#dce5f9'),axis.text=element_text(color='#dce5f9'),panel.grid=element_line(color='#2b3853'),plot.background=element_rect(fill='#121b32',color=NA),panel.background=element_rect(fill='#121b32',color=NA))
  },bg='#121b32')
  output$error_caption <- renderText({r <- result();switch(input$error_region,alpha=sprintf('Type I error: %.1f%% of null-distribution area is in the rejection region.',100*input$alpha),correct=sprintf('Correct non-rejection under H₀: %.1f%%.',100*(1-input$alpha)),beta=sprintf('Type II error at the specified effect: %.1f%%.',100*(1-r$power)),power=sprintf('Power at the specified effect: %.1f%%.',100*r$power))})
  output$comparison <- renderTable({a <- result(); b <- calc_power(input$compare_n,input$d,input$alpha,input$test,input$alternative);data.frame(Measure=c('n (per group for two-sample)','Total N','Statistical power'),Current=c(input$n,a$total,sprintf('%.1f%%',100*a$power)),Compared=c(input$compare_n,b$total,sprintf('%.1f%%',100*b$power)))},striped=TRUE)
  challenge_started <- reactiveVal(FALSE)
  observeEvent(input$challenge,{challenge_started(TRUE);showModal(modalDialog(title='Challenge: reach 80% power',p('Change sample size until the current design reaches at least 80% power. Keep the assumed effect and alpha fixed.'),easyClose=TRUE,footer=modalButton('Let me try')))})
  output$challenge_text <- renderText({if(!challenge_started()) return('');if(result()$power>=.8) 'Challenge completed: the current configuration has at least 80% power!' else 'Challenge in progress: increase the sample size to reach at least 80% power.'})
  sensitivity_d <- reactive({
    f <- function(d) calc_power(input$n,d,input$alpha,input$test,input$alternative)$power-input$target
    direction <- if(input$alternative=='less') -1 else 1
    if(f(0)>=0) return(0)
    if(f(direction*20)<0) return(NA_real_)
    uniroot(function(z) f(direction*z),c(0,20),tol=1e-7)$root*direction
  })
  output$sensitivity_result <- renderText({d <- sensitivity_d();if(is.na(d)) 'Target not reached for effect magnitudes up to 20.' else sprintf('Minimum effect for %.0f%% target power in the selected direction: Cohen’s d = %.3f',100*input$target,d)})
  output$download_report <- downloadHandler(filename=function()'PowerLab_analysis.txt',content=function(file){
    r <- result(); writeLines(c('PowerLab — prospective power analysis',paste('Date:',as.character(Sys.Date())),paste('Test:',if(input$test=='one')'One-sample t-test' else 'Independent two-sample pooled-variance t-test (equal n and variances)'),paste('Alternative:',input$alternative),paste('n:',input$n),paste('Total N:',r$total),paste('Assumed signed Cohen d:',input$d),paste('Alpha:',input$alpha),paste('Power:',sprintf('%.6f',r$power)),paste('Target power:',input$target),paste('Suggested n:',if(is.na(recommendation()))'Not reached up to 5000' else recommendation()),paste('df:',r$df),paste('Noncentrality parameter:',r$delta),'Assumptions: normal independent observations, fixed true effect size; pooled equal-variance and equal n for two groups.'),file)
  })
  output$download_csv <- downloadHandler(filename=function()'PowerLab_power_curve.csv',content=function(file)write.csv(curve_data(),file,row.names=FALSE))
  observeEvent(input$curve_click,{
    click <- input$curve_click
    if(is.null(click$x)||is.null(click$y)) return()
    if(abs(click$x-input$n)<5 && abs(click$y-result()$power)<.08) show_explanation('current')
    else if(abs(click$y-input$target)<.035) show_explanation('targetline')
    else {
      closest_n <- max(3,min(200,round(click$x)))
      expected <- calc_power(closest_n,input$d,input$alpha,input$test,input$alternative)$power
      if(abs(click$y-expected)<.09) show_explanation('curve',sprintf('At n = %d, predicted power = %.1f%%.',closest_n,100*expected))
    }
  })
}
shinyApp(ui,server)
