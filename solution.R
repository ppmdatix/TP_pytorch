%install.packages("bnlearn")
%source("http://bioconductor.org/biocLite.R")
%biocLite(c("graph", "Rgraphviz"))

library("bnlearn")

# charge la base de données alarm
data("alarm")

# infos sur les colonnes (nos variables)
ncol(alarm)
colnames(alarm)

# infos sur les lignes (nos observations)
nrow(alarm)
rownames(alarm)

# dix premières lignes, 5 premières colonnes
alarm[1:10, 1:5]

# lignes 3, 5, 1, colonnes "ANES", "HIST" et "MINV"
alarm[c(3, 5, 1), c("ANES", "HIST", "MINV")]

# 1. Qu'obtient-on avec les commandes suivantes?

table(alarm[, "PAP"])
plot(alarm[, "PAP"])
prop.table(table(alarm[, "PAP"]))

table(alarm[, "SHNT"])
plot(alarm[, "SHNT"])
prop.table(table(alarm[, "SHNT"]))

ct = table(alarm[, c("PAP", "SHNT")])
prop.table(ct)
prop.table(ct, margin = 1)
prop.table(ct, margin = 2)

res = ci.test(x = "PAP", y = "SHNT", z = as.character(NULL), data = alarm, test = "mi")
res$statistic
res$p.value
res

res = ci.test(x = "PAP", y = "SHNT", z = "PMB", data = alarm, test = "mi")
res$statistic
res$p.value
res

# 2. Relations à vérifier

ci.test(x = "STKV", y = "HR", data = alarm, test = "mi")
ci.test(x = "STKV", y = "HR", z = "CO", data = alarm, test = "mi")
ci.test(x = "HR", y = "CO", data = alarm, test = "mi")
ci.test(x = "HR", y = "CO", z = "STKV", data = alarm, test = "mi")
ci.test(x = "CO", y = "STKV", data = alarm, test = "mi")
ci.test(x = "CO", y = "STKV", z = "HR", data = alarm, test = "mi")

# 3. Inspection de la relation entre STKV et HR sachant CO

mask = rep(TRUE, nrow(alarm))
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x)")

mask = alarm[, "CO"] == "HIGH"
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x,z=HIGH)")

mask = alarm[, "CO"] == "LOW"
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x,z=LOW)")

mask = alarm[, "CO"] == "NORMAL"
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x,z=NORMAL)")

# 4. Inférence

# apprentissage de structure
bn = hc(alarm)
graphviz.plot(bn)

# apprentissage de parametres
bn = bn.fit(bn, data = alarm, method = "bayes")
bn[["CO"]]

# 5. Requetes approchées
cpquery(bn, event = (STKV == "HIGH"), evidence = (HR == "LOW"))
cpquery(bn, event = (STKV == "HIGH"), evidence = (HR == "LOW" & CO == "LOW"))

# 6. Requetes exactes
source("includes.R")

p = exact.dist(bn, event = c("STKV", "HR", "CO"), evidence = TRUE)
sum(p["HIGH", "LOW", ]) / sum(p[, "LOW", ])
sum(p["HIGH", "LOW", "LOW"]) / sum(p[, "LOW", "LOW"])

p = exact.dist(bn, event = c("INT", "APL"), evidence = TRUE)

# 6. Do-calculus

# p(y|x) et p(y|do(x)) avec y=HYP et x=STKV
p = exact.dist(bn, event = c("HYP", "STKV"), evidence = TRUE)
p.y.x = prop.table(p, margin = 2)

p = exact.dist(bn, event = c("HYP", "STKV", "LVV"), evidence = TRUE)
p.z = margin.table(p, margin = 3)
p.y.xz = prop.table(p, margin = c(2, 3))
p.y.do.x = margin.table(p.y.xz * rep(p.z, each=prod(dim(p.y.xz)[-3])), margin = c(1, 2))

p.y.x
p.y.do.x

# 9. Algorithme PC

vars = colnames(alarm)

# graphe complet
g = empty.graph(vars)
for (x in vars) {
  for (y in setdiff(vars, x)) {
    g = set.edge(g, from = x, to = y)
  }
}
graphviz.plot(g)

# squelette
alpha = 0.001
z.xy = sapply(vars, function(v) {list()})
for(m in 0:length(vars)) {
  
  x.done = NULL
  for (x in vars) {
    
    x.done = c(x.done, x)
    x.nbrs = g$nodes[[x]]$nbr
    
    for (y in setdiff(x.nbrs, x.done)) {
      y.nbrs = g$nodes[[y]]$nbr
      
      z.cands = array(NA, dim = c(m, 0))
      if(length(x.nbrs) > m) {
        z.cands = cbind(z.cands, combn(setdiff(x.nbrs, y), m))
      }
      if(length(y.nbrs) > m) {
        z.cands = cbind(z.cands, combn(setdiff(y.nbrs, x), m))
      }
      
      for (i in (0:ncol(z.cands))[-1]) {
        z = z.cands[, i]
        
        cat("testing", x, "indep", y, "given", paste(z, collapse = ", "), "\n")
        res = ci.test(x = x, y = y, z = z, data = alarm, test = "mi")
        
        if(res$p.value > alpha) {
          cat("dropping edge\n")
          g = drop.edge(g, from = x, to = y)
          x.nbrs = setdiff(x.nbrs, y)
          z.xy[[x]][[y]] = z
          break
        }
      }
    }
  }
}
skeleton = g

graphviz.plot(g)

# graphe (presque) dirigé

g = skeleton
x.done = NULL
for (x in vars) {
  
  x.done = c(x.done, x)
  x.nbrs = g$nodes[[x]]$nbr
  
  for (w in setdiff(x.nbrs, x.done)) {
    w.nbrs = g$nodes[[w]]$nbr
    
    for (y in setdiff(w.nbrs, c(x.done, x.nbrs))) {
      
      if (!(w %in% z.xy[[x]][[y]])) {
        if (w %in% g$nodes[[x]]$parents) {
          cat("warning, inconsistent arc", x, "<->", w, "\n")
        }
        if (w %in% g$nodes[[y]]$parents) {
          cat("warning, inconsistent arc", y, "<->", w, "\n")
        }
        
        g = set.arc(g, from = x, to = w)
        g = set.arc(g, from = y, to = w)
      }
    }
  }
}
graphviz.plot(g)
