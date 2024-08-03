rm(list = ls())

dir <- ""

is <- 1:9
datas <- c("prostate", "glucose")
cols <-  hcl.colors(5, palette = "Zissou 1", alpha = 0.7)[c(1, 2, 5)]

stim <- sdiv <- sess <- array(NA, dim = c(length(datas), length(is), 10, 3))
supl <- array(NA, dim = c(length(datas), length(is), 1000, 10, 3))

qv <- array(NA, dim = c(length(is)))

for (data in datas) {
  for (i in is) {
    load(paste(dir, "Out/Data/simple_", data, "_q", i, ".RData", sep = ""))
    sess[which(data == datas), which(i == is), , ] <- ess[ , c(2, 1, 3)]
    stim[which(data == datas), which(i == is), , ] <- tim[ , c(2, 1, 3)]
    sdiv[which(data == datas), which(i == is), , ] <- div[ , c(2, 1, 3)]
    supl[which(data == datas), which(i == is), , , 1] <- upl2
    supl[which(data == datas), which(i == is), , , 2] <- upl1
    supl[which(data == datas), which(i == is), , , 3] <- upl3
    qv[which(i == is)] <- q
  }
}

mupl <- apply(supl, c(1, 2, 4, 5), mean)

pdf(paste(dir, "Out/Figures/res.pdf", sep = ""), width = 6, height = 4,
    family = "Times")
lo <- cbind(c(1, 2), c(3, 4), c(5, 6))
layout(lo)
par(mar = c(1, 4, 1, 0))
par(oma = c(3, 2, 3, 1))

plot.datas <- c("prostate", "glucose")
for (data in plot.datas) {
  sub <- t(apply(aperm(mupl[which(data == datas), , , ], c(3, 2, 1)), 2, "c"))
  plot(1:ncol(sub), sub[1, ],
       col = rep(cols, length(qv)), ylim = range(sub) + c(0, 1),
       axes = FALSE, type = "n", xlab = "", ylab = "", main = "")
  box()
  if (data == plot.datas[1]) {
    axis(2, c(-90, -60))
    mtext("Prostate", 2, line = 4.5, cex = 0.8)
  } else {
    axis(2, c(-400, -200))
    mtext("Glucose", 2, line = 4.5, cex = 0.8)
  }
  for (j in 1:nrow(sub)) {
    points(1:ncol(sub), sub[j, ], col = rep(cols, length(qv)), pch = 16)
  }
  abline(v = seq(0.5, ncol(sub) + 3, by = 3), lty = 3, col = "gray")
  if (data == plot.datas[1]) {
    mtext(expression(widehat(E)*"["*frac(1,2*sigma^2)*"||y-X"*z[2]*"||"[2]^2+lambda*"||"*z[2]*"||"[q]^q*"| y]"), 3)
    legend("bottomright", col = cols[c(3, 1, 2)], pch = 16, legend = c("Naive (1)", "Centered (2)", "Non-Centered (3)"),
           bg = "white")
    axis(1, seq(2, ncol(sub), by = 3), qv)
  } else {
    axis(1, seq(2, ncol(sub), by = 3), qv)
    mtext("q", 1, line = 2)
  }
  mtext(expression(widehat(E)*"["*frac(1,2*sigma^2)*"||y-X"*z[2]*"||"[2]^2+lambda*"||"*z[2]*"||"[q]^q*"| y]"), 2, line = 2, cex = 0.8)
}

for (data in plot.datas) {
  sub <- t(apply(aperm(sess[which(data == datas), , , ], c(3, 2, 1)), 2, "c"))
  plot(1:ncol(sub), sub[1, ],
       col = rep(cols, length(qv)), ylim = range(sub) + c(0, 1),
       axes = FALSE, type = "n", xlab = "", ylab = "", main = "")
  box()
  for (j in 1:nrow(sub)) {
    points(1:ncol(sub), sub[j, ], col = rep(cols, length(qv)), pch = 16)
  }
  if (data == plot.datas[1]) {
    mtext("Min. ESS", 3, line = 1)
    axis(1, seq(2, ncol(sub), by = 3), qv)
  } else {
    axis(1, seq(2, ncol(sub), by = 3), qv)
  }
  abline(v = seq(0.5, ncol(sub), by = 3), lty = 3, col = "gray")
  if (data == plot.datas[1]) {
    axis(2, c(200, 800))
  } else {
    axis(2, c(200, 800))
    mtext("q", 1, line = 2)
  }
  mtext("Min. ESS", 2, line = 2, cex = 0.8)
}

for (data in plot.datas) {
  sub <- log(t(apply(aperm(stim[which(data == datas), , , ], c(3, 2, 1)), 2, "c")))
  plot(1:ncol(sub), sub[1, ],
       col = rep(cols, length(qv)), ylim = range(sub) + c(0, 1),
       axes = FALSE, type = "n", xlab = "", ylab = "", main = "")
  box()
  for (j in 1:nrow(sub)) {
    points(1:ncol(sub), sub[j, ], col = rep(cols, length(qv)), pch = 16)
  }
  if (data == plot.datas[1]) {
    mtext("log(Time)", 3, line = 1)
    axis(1, seq(2, ncol(sub), by = 3), qv)
  } else {
    axis(1, seq(2, ncol(sub), by = 3), qv)
  }
  abline(v = seq(0.5, ncol(sub), by = 3), lty = 3, col = "gray")
  if (data == plot.datas[1]) {
    axis(2, c(-1, 2))
  } else {
    axis(2, c(0, 3))
    mtext("q", 1, line = 2)
  }
  mtext("log(Time)", 2, line = 2, cex = 0.8)
}
dev.off()

pdf(paste(dir, "Out/Figures/div.pdf", sep = ""), width = 6, height = 2,
    family = "Times")
lo <- rbind(c(1, 2))
layout(lo)
par(mar = c(1, 2, 1, 0))
par(oma = c(2, 2, 0, 1))
for (data in plot.datas) {
  sub <- t(apply(aperm(sdiv[which(data == datas), , , ], c(3, 2, 1)), 2, "c"))
  plot(1:ncol(sub), sub[1, ], ylim = range(sub) + c(0, 1),
       col = rep(cols, length(qv)),
       axes = FALSE, type = "n", xlab = "", ylab = "", main = "")
  box()
  for (j in 1:nrow(sub)) {
    points(1:ncol(sub), sub[j, ], col = rep(cols, length(qv)), pch = 16)
  }
  if (data == plot.datas[1]) {
    mtext("Prostate", 3, line = 0)
    mtext("Divergences", 2, line = 3)
    mtext("after Warmup", 2, line = 2)
    axis(1, seq(2, ncol(sub), by = 3), qv)
  } else {
    mtext("Glucose", 3, line = 0)
    axis(1, seq(2, ncol(sub), by = 3), qv)
  }
  abline(v = seq(0.5, ncol(sub), by = 3), lty = 3, col = "gray")
  if (data == plot.datas[1]) {
    axis(2, c(0, 300, 600))
    legend("topright", col = cols[c(3, 1, 2)], pch = 16, legend = c("Naive (1)", "Centered (2)", "Non-Centered (3)"),
           bg = "white", cex = 0.8)
  } else {
    axis(2, c(0, 400, 800))
  }
  mtext("q", 1, line = 1.9)
}
dev.off()


cols <-  hcl.colors(5, palette = "Zissou 1", alpha = 1)[c(1, 2, 5)]

pdf(paste(dir, "Out/Figures/kdens.pdf", sep = ""), width = 6, height = 4,
    family = "Times")
lo <- rbind(c(1, 4), c(2, 5), c(3, 6))
layout(lo)
par(mar = c(1, 2, 3, 0))
par(oma = c(5, 3, 0, 1))
for (data in plot.datas) {

  dupl1 <- lapply(1:10, function(chain, data) {
    density(supl[which(data == datas), 1, , chain, 1])
  }, data = data)
  dupl2 <- lapply(1:10, function(chain, data) {
    density(supl[which(data == datas), 1, , chain, 2])
  }, data = data)
  dupl3 <- lapply(1:10, function(chain, data) {
    density(supl[which(data == datas), 1, , chain, 3])
  }, data = data)

  for (mod in c(3, 1, 2)) {
plot(dupl1[[1]]$x, dupl1[[1]]$y,
     ylim = c(0, 0.4),
     xlim = range(unlist(lapply(dupl1, function(x) {x$x})),
                  unlist(lapply(dupl2, function(x) {x$x})),
                  unlist(lapply(dupl3, function(x) {x$x}))),
     type = "n", axes = FALSE, xlab = "", ylab = "")
    mtext(ifelse(mod == 1, "Centered (2)",
                 ifelse(mod == 2, "Non-Centered (3)", "Naive (1)")), 3,
          line = 0.5)
for (chain in 1:10) {
  if (mod == 1) {
    lines(dupl1[[chain]]$x, dupl1[[chain]]$y, col = cols[1])
  } else if (mod == 2) {
    lines(dupl2[[chain]]$x, dupl2[[chain]]$y, col = cols[2])
  } else {
    lines(dupl3[[chain]]$x, dupl3[[chain]]$y, col = cols[3])
  }
}
if (data == plot.datas[1]) {
  axis(1, c(-200, -110, -90, 0), line = 0)
  axis(2, c(-1, 0.1, 0.3, 1))
  mtext("Kernel Density", 2, line = 3)
  mtext("Estimate", 2, line = 2)
} else {
  axis(1, c(-600, -500, -200, 0), line = 0)
  axis(2, c(-1, 0.1, 0.3, 1), rep("", 4))
}
    if (mod == 2) {
      mtext(expression(frac(1,2*sigma^2)*"||y-X"*z[2]*"||"[2]^2+lambda*"||"*z[2]*"||"[q]^q), 1,
            line = 4.5)
    }
    if (mod == 3) {
      mtext(ifelse(data == "prostate", "Prostate", "Glucose"), 3,
            line = 1.7)
    }
    axis(3, c(-1000, 1000))
    axis(4, c(-1, 1))
  }

}
dev.off()
