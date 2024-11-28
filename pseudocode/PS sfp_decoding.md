# Pseudocode for behavioral analyses

##SVM and basic decoding
1. **Initialize Parameters and Load Data**
   - Set the number of odors and data window size.
   - Define directories containing sniff feature data.
   - Load behavioral data needed for analysis.

2. **Correlation Analysis (Optional)**
   - If enabled:
     - For each subject:
       - Load sniff features and trial onset times.
       - Organize trials by odor identity.
       - Compute correlation matrices of sniff features.
       - Compare correlations for trials of the same odor versus different odors.
     - Visualize results with a bar plot.

3. **Support Vector Machine (SVM) Decoding**
   - For each subject:
     - Load sniff features and trial onset times.
     - Organize trials by odor identity.
     - Preprocess sniff features (e.g., select relevant features, handle missing data).
     - Perform Principal Component Analysis (PCA) to reduce dimensionality.
     - Train an SVM classifier to predict odor identity from sniff features.
     - Evaluate classifier performance (accuracy).
     - Analyze behavioral correlations between actual and predicted odors for misclassified trials.
     - Perform statistical tests to assess significance.
   - Visualize SVM performance and behavioral correlations.
   - Save results for future analysis.

