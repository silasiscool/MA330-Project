# RUN THIS FILE
# REQUIRED LIBRARIES: matplotlib, sympy, numpy, scipy

# Modules
import matplotlib.pyplot as plt
from sympy import pi, N, cos, sin, Piecewise, lambdify
import numpy as np
import fourier
from errorCalc import NRMSE

# Paremeters
errorThresholds = [0.10, 0.05, 0.01]
P = 2*pi # Period
A = 1 # Amplitude

# Vars
t = fourier.t
n = fourier.n
L = P/2
tSpace = np.linspace(-N(P)/2, N(P)/2, 100) # time values to plot over

# Waveform functions (with amplitude A and period P)
saw = lambda t : 2*A/P * t # Sawtooth wave defined over [-P/2, P/2]
saw_sym = saw(t)
square_sym = Piecewise((A, t >= 0), (-A, True))
square = lambdify(t, square_sym)








# Sawtooth Transform
a0 = fourier.a0(saw_sym, P)
an = fourier.an(saw_sym, P)
bn = fourier.bn(saw_sym, P)

print(f'Sawtooth: a0: {a0}, an: {an}, bn: {bn}')

sawTerms = [a0]
errors = []
i = 1
thresholdNum = 0
thresholdNs = []
thresholdSeries = []
while True:
    anTerm = an.subs(n,i) * cos(i*pi*t/L)
    bnTerm = bn.subs(n,i) * sin(i*pi*t/L)
    sawTerms.append(anTerm)
    sawTerms.append(bnTerm)
    
    sawSeries = sum(sawTerms)
    error = NRMSE(saw_sym, sawSeries, t, P, A)
    errors.append(error)

    if error <= errorThresholds[thresholdNum]:
        thresholdNs.append(i)
        thresholdSeries.append(sawSeries)
        thresholdNum+=1
        if len(errorThresholds) <= thresholdNum:
            break
    i+=1


# Create figure for sawtooth
fig = plt.figure(1)
fig.suptitle("Sawtooth Waveform Fouirer Transform")

# Plot original function and each series
for i in range(len(thresholdSeries)):
    ax = plt.subplot(2,len(thresholdSeries),i+1)
    waveformPlt = plt.plot(tSpace, [saw_sym.subs(t, tval) for tval in tSpace])
    seriesPlt = plt.plot(tSpace, [thresholdSeries[i].subs(t, tval) for tval in tSpace])

    ax.set_title(f'NRMSE threshold: {errorThresholds[i]}, {thresholdNs[i]} Terms')

# Plot error vs # terms
ax = plt.subplot(2,len(thresholdSeries),len(thresholdSeries)+1)
plt.plot(range(1, len(errors)+1), errors)

ax.set_title('NRMSE vs # Terms')
ax.set_ylabel('NRSME')
ax.set_xlabel('# terms')

plt.tight_layout()







# Square Wave Transform
a0 = fourier.a0(square_sym, P)
an = fourier.an(square_sym, P)
bn = fourier.bn(square_sym, P)

print(f'Square: a0: {a0}, an: {an}, bn: {bn}')

squareTerms = [a0]
errors = []
i = 1
thresholdNum = 0
thresholdNs = []
thresholdSeries = []
while True:
    anTerm = an.subs(n,i) * cos(i*pi*t/L)
    bnTerm = bn.subs(n,i) * sin(i*pi*t/L)
    squareTerms.append(anTerm)
    squareTerms.append(bnTerm)
    
    squareSeries = sum(squareTerms)
    error = NRMSE(square_sym, squareSeries, t, P, A)
    errors.append(error)

    if error <= errorThresholds[thresholdNum]:
        thresholdNs.append(i)
        thresholdSeries.append(squareSeries)
        thresholdNum+=1
        if len(errorThresholds) <= thresholdNum:
            break
    i+=1


# Create figure for square
fig = plt.figure(2)
fig.suptitle("Square Waveform Fouirer Transform")

# Plot original function and each series
for i in range(len(thresholdSeries)):
    ax = plt.subplot(2,len(thresholdSeries),i+1)
    plt.plot(tSpace, [square_sym.subs(t, tval) for tval in tSpace])
    plt.plot(tSpace, [thresholdSeries[i].subs(t, tval) for tval in tSpace])

    ax.set_title(f'NRMSE threshold: {errorThresholds[i]}, {thresholdNs[i]} Terms')

# Plot error vs # terms
ax = plt.subplot(2,len(thresholdSeries),len(thresholdSeries)+1)
plt.plot(range(1, len(errors)+1), errors)

ax.set_title('NRMSE vs # Terms')
ax.set_ylabel('NRSME')
ax.set_xlabel('# terms')

plt.tight_layout()








# Show plots
plt.show()
