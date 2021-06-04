low = 75;
high = 200;

xrange = tlen[low:high];
yrange = ylen[shock_locations[low:high]];

xrange = @. u"ns"(xrange); xrangeu = ustrip(xrange);
yrange = @. u"μm"(yrange); yrangeu = ustrip(yrange);

# Compute slope
data = DataFrame(y = yrangeu, x = xrangeu);
slope = lm(@formula(y ~ x), data);

b, m = coef(slope);

f(x) = @. (m*x + b) * u"μm";

slp = m * u"μm/ns";

pl = plot(xrange, [yrange, f.(xrangeu)], labels = ["Grad vals" "Regression"], title = "Slope = $slp");
png(pl, shotstring*"slp")