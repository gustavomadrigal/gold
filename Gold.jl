using LaTeXStrings
using Plots
using DelimitedFiles
using Statistics 

function lennard_jones_potential(r, ε, σ)
    return 4ε * ((σ / r)^12 - (σ / r)^6)
end

# Define parameters
ε = 0.584749  # Depth of the potential energy well
σ = 2.593  # Finite distance at which the potential is zero

# Generate distance values
r_values = range(1, stop=10.0, length=10000)

# Calculate potential energy values
V_values = lennard_jones_potential.(r_values, ε, σ)

# Plot the Lennard-Jones potential
plot(r_values, V_values, xlabel=L"r(\AA)", ylabel="V (eV)", label="Lennard-Jones Potential", legend=true)

ylims!(-6.6, 10)
xlims!(1.3, 9)

#hline!([0], linestyle=:dash, color=:black, label="Zero Energy Line")
hline!([-0.584749], linestyle=:dash, color=:red, label="Minimum Energy")

#savefig("figure.pdf")


data = readdlm("rdf_1000K.txt")

# Extract x and y values
x = data[:, 2]
y = data[:, 3]

# Create a plot
plot(x, y, xlabel =L"r (\AA)", ylabel = " g(r)", title = "Radial Distribution Function",label = "1000 K")
ylims!(0,50)
savefig("rdf_1000.pdf")

data = readdlm("rdf_1500K.txt")

# Extract x and y values
x = data[:, 2]
y = data[:, 3]

# Create a plot
plot(x, y, xlabel =L"r (\AA)", ylabel = " g(r)", title = "Radial Distribution Function",label = "1500 K")
ylims!(0,50)
savefig("rdf_1500.pdf")

data = readdlm("rdf_2000K.txt")

# Extract x and y values
x = data[:, 2]
y = data[:, 3]

# Create a plot
plot(x, y, xlabel =L"r (\AA)", ylabel = " g(r)", title = "Radial Distribution Function",label = "2000 K")
ylims!(0,50)
savefig("rdf_2000.pdf")

data = readdlm("rdf_2500K.txt")

# Extract x and y values
x = data[:, 2]
y = data[:, 3]

# Create a plot
plot(x, y, xlabel =L"r (\AA)", ylabel = " g(r)", title = "Radial Distribution Function",label = "2500 K")
ylims!(0,50)
savefig("rdf_2500.pdf")

data = readdlm("rdf_3000K.txt")

# Extract x and y values
x = data[:, 2]
y = data[:, 3]

# Create a plot
plot(x, y, xlabel =L"r (\AA)", ylabel = " g(r)", title = "Radial Distribution Function",label = "3000 K")
ylims!(0,50)
savefig("rdf_3000.pdf")

data1 = readdlm("results_500.txt")
data2 = readdlm("results_1000.txt")
data3 = readdlm("results_1500.txt")

# Extract x and y values
x1 = data1[:, 1]
y1 = data1[:, 11]*(1/6)

x2 = data2[:, 1]
y2 = data2[:, 11]*(1/6)

x3 = data3[:, 1]
y3 = data3[:, 11]*(1/6)

# Create a plot
plot(x1, y1, xlabel ="t (ps)", ylabel =L" \textrm{MSD (\AA^{2})}",label = "500 K")
plot!(x2,y2,label = "1000 K")
plot!(x3,y3,label = "1500 K")
#xlims!(0,1100)
#ylims!(0,0.025)
#savefig("msd.pdf")

data11 = readdlm("results_500_linear.txt")

# Extract x and y values
x11 = data11[:, 1]
y11 = (data11[:, 11])*(10^-10/6)

plot(x11, y11, xlabel ="t (ps)",  ylabel =L" \textrm{MSD (m^{2})}",label = "500 K")

# Insert text at a specific point on the plot
annotate!(750, 4*10^-13, text("MSD = 1.47e-15*t - 3.65e-13", 10, :left, :center, :black))

# Perform linear regression
n = length(x11)
X = [ones(n) x11]
coeffs = X \ y11  # Solve the normal equations
intercept = coeffs[1]
slope = coeffs[2]

# Generate the fitted line
x_fit = range(minimum(x11), maximum(x11), length = 10)
y_fit = slope * x_fit .+ intercept

# Plot the fitted line
scatter!(linestyle=:dash,x_fit, y_fit, label = "Linear Fit")
#savefig("msd_500.pdf")

# Calculate R-squared
y_mean = mean(y1)
y_pred = slope * x11 .+ intercept
SS_total = sum((y11 .- y_mean).^2)
SS_residual = sum((y11 .- y_pred).^2)
R_squared = 1 - SS_residual / SS_total

# Calculate fitting equation
fit_equation = "y = $(round(slope; digits=40))*x + $(round(intercept; digits=40))"

# Format R-squared value
r_squared_text = "R² = $(round(R_squared; digits=10))"

println("Fitting Equation: $fit_equation")
println("R-squared: $r_squared_text")

data22 = readdlm("results_1000_linear.txt")

# Extract x and y values
x22 = data22[:, 1]
y22 = data22[:, 11]*(10^-10/6)

plot(x22, y22, xlabel ="t (ps)", ylabel =L" \textrm{MSD (m^{2})}",label = "1000 K")

# Insert text at a specific point on the plot
annotate!(750, 5*10^-13, text("MSD = 1.84e-15*t - 4.54e-13", 10, :left, :center, :black))

# Perform linear regression
n = length(x11)
X = [ones(n) x11]
coeffs = X \ y22  # Solve the normal equations
intercept = coeffs[1]
slope = coeffs[2]

# Generate the fitted line
x_fit = range(minimum(x22), maximum(x22), length = 10)
y_fit = slope * x_fit .+ intercept

# Plot the fitted line
scatter!(linestyle=:dash,x_fit, y_fit, label = "Linear Fit")
#savefig("msd_1000.pdf")

# Calculate R-squared
y_mean = mean(y22)
y_pred = slope * x22 .+ intercept
SS_total = sum((y22 .- y_mean).^2)
SS_residual = sum((y22 .- y_pred).^2)
R_squared = 1 - SS_residual / SS_total

# Calculate fitting equation
fit_equation = "y = $(round(slope; digits=40))*x + $(round(intercept; digits=40))"

# Format R-squared value
r_squared_text = "R² = $(round(R_squared; digits=10))"

println("Fitting Equation: $fit_equation")
println("R-squared: $r_squared_text")

data33 = readdlm("results_1500_linear.txt")

# Extract x and y values
x33 = data33[:, 1]
y33 = data33[:, 11]*(10^-10/6)

plot(x33, y33, xlabel ="t (ps)", ylabel =L" \textrm{MSD (m^{2})}",label = "1500 K")

# Insert text at a specific point on the plot
annotate!(750,5* 10^-13, text("MSD = 2.74e-15*t - 6.75e-13", 10, :left, :center, :black))

# Perform linear regression
n = length(x33)
X = [ones(n) x33]
coeffs = X \ y33  # Solve the normal equations
intercept = coeffs[1]
slope = coeffs[2]

# Generate the fitted line
x_fit = range(minimum(x22), maximum(x22), length = 10)
y_fit = slope * x_fit .+ intercept

# Plot the fitted line
scatter!(linestyle=:dash,x_fit, y_fit, label = "Linear Fit")
#savefig("msd_1500.pdf")

# Calculate R-squared
y_mean = mean(y33)
y_pred = slope * x33 .+ intercept
SS_total = sum((y33 .- y_mean).^2)
SS_residual = sum((y33 .- y_pred).^2)
R_squared = 1 - SS_residual / SS_total

# Calculate fitting equation
fit_equation = "y = $(round(slope; digits=40))*x + $(round(intercept; digits=40))"

# Format R-squared value
r_squared_text = "R² = $(round(R_squared; digits=10))"

println("Fitting Equation: $fit_equation")
println("R-squared: $r_squared_text")

d = [1.0e-15,2.0e-15,3.0e-15]
plot(1 ./T,log10.(d), xlabel=L"\textrm{1/T (K^{-1})}", ylabel =L"\textrm{Log[D] (m^{2}/s)}", label = "Simulation Data")

# Insert text at a specific point on the plot
annotate!(0.0011, -14.9, text("Log(D) = -344.73*(1/T) + -14.32", 10, :left, :center, :black))

# Perform linear regression
n = length(1 ./T)
X = [ones(n) 1 ./T]
coeffs = X \ log10.(d)  # Solve the normal equations
intercept = coeffs[1]
slope = coeffs[2]

# Generate the fitted line
x_fit = range(minimum(1 ./T), maximum(1 ./T), length = 10)
y_fit = slope * x_fit .+ intercept

# Plot the fitted line
scatter!(linestyle=:dash,x_fit, y_fit, label = "Linear Fit")
#savefig("ea.pdf")

# Calculate R-squared
y_mean = mean(d)
y_pred = slope * (1 ./T) .+ intercept
SS_total = sum((log10.(d) .- y_mean).^2)
SS_residual = sum((log10.(d) .- y_pred).^2)
R_squared = 1 - SS_residual / SS_total

# Calculate fitting equation
fit_equation = "y = $(round(slope; digits=40))*x + $(round(intercept; digits=40))"

# Format R-squared value
r_squared_text = "R² = $(round(R_squared; digits=10))"

println("Fitting Equation: $fit_equation")
println("R-squared: $r_squared_text")


