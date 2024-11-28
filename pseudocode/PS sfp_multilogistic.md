Pseudocode for multilogistic model

1. **Initialize Parameters and Settings**
   - Set flags and parameters (e.g., `grp_`, `num_descrip`, `settings_.pcamaker`).
   - Define file paths and directories.

2. **Load Sniff Labels**
   - Load sniff feature labels from 'snifflabels.mat'.
   - Adjust the label list as needed.

3. **Load Behavioral Data**
   - If grouping (`grp_`) is enabled:
     - Set the number of odors and sample window size.
     - Define directories for each subject's data.
     - Load behavioral ratings from 'NEMO_perceptual2.mat'.
     - Initialize variables to store results.

4. **Process Data for Each Subject**
   FOR each subject in the data directories:
   a. **Load Sniff Features**
      - Load sniff feature matrices for the subject.
      - Concatenate feature matrices into a single matrix.
      - Prune and preprocess the feature matrix (e.g., select relevant features, handle NaNs, normalize).
   
   b. **Load Trial Onset Information**
      - Load onset times for the subject.
      - Create grouping vectors and matrices to organize trial data.
   
   c. **Optional: Perform PCA**
      - If PCA is enabled, reduce dimensionality of the feature matrix.
   
   d. **Multilogistic Regression for Each Perceptual Descriptor**
      FOR each perceptual descriptor:
      - Retrieve behavioral ratings for the descriptor.
      - Discretize ratings into bins for classification.
      - Fit a multilogistic regression model using cross-validation.
      - Store the effectiveness and p-values of the model.
      - If PCA was performed, adjust weights accordingly.
   
   e. **Visualize Results**
      - Generate a heatmap or image to visualize regression weights for the subject.

5. **Aggregate and Analyze Results Across Subjects**
   - Adjust perceptual descriptor labels to align across subjects.
   - Compute mean effectiveness scores and sort descriptors accordingly.
   - Plot aggregated results with error bars to display variability across subjects.
   - Add appropriate labels and formatting to the plots.

6. **Optional: Perform FDR Correction**
   - Apply False Discovery Rate correction to p-values if needed.

