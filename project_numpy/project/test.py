import numpy as np
import matplotlib.pyplot as plt

value = int(input('Enter a number: '))

# Generate random data
data = np.random.randn(value)

# Create histogram
plt.hist(data, bins=20)

# Set labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram')

# Display the histogram
plt.show()