# matplotlib_Pymaceuticals
Real-world data analysis on the fictional company Pymaceuticals using pandas and matplotlib.
![Laboratory](Images/Laboratory.jpg)

The tasks, python code, results and figures are shown below. Here are three observations made from observing the data and plots:
* An obvious observation is the breakdown of male versus female subjects. The analysis is almost 50/50. Having an even gender distribution is a good thing for testing drugs as a drug may have different effects on different genders.

* Looking at the box plot, it is clear that Capomulin and Ramicane are the leading regimen options. The final tumor size ranges for these are smaller than that of Infubinol and Ceftamin.

* From running a linear regression on mouse weight versus average tumor volume, we find the r-value is 0.84. An r-value close to one means there is a good correlation between these two parameters. The correlation is seen visually by observing how close the regression line is to the scatter plot points.

# Tasks, Code, Results and Figures
## Imports, Load and Clean Data
```python
# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import linregress

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata_df = pd.read_csv(mouse_metadata_path)
study_results_df = pd.read_csv(study_results_path)

# Combine the data into a single dataset
combine_data_df = pd.merge(mouse_metadata_df, study_results_df, how="inner")

# Display the data table for preview
combine_data_df.head()
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Drug Regimen</th>
      <th>Sex</th>
      <th>Age_months</th>
      <th>Weight (g)</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>0</td>
      <td>45.000000</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>5</td>
      <td>38.825898</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>10</td>
      <td>35.014271</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>15</td>
      <td>34.223992</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>20</td>
      <td>32.997729</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>

```python
# Checking the number of mice.
unique_mice_id_count = len(combine_data_df["Mouse ID"].unique())
print(f"Unique Mice Id's: {unique_mice_id_count}")
```
Unique Mice Id's: 249

```python
# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mice_df = combine_data_df[combine_data_df.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"]
duplicate_mice_df.head()
```
909    g989  
911    g989  
913    g989  
915    g989  
917    g989  
Name: Mouse ID, dtype: object
```python
# Optional: Get all the data for the duplicate mouse ID. 
duplicate_mice_df = combine_data_df[combine_data_df.duplicated(subset=["Mouse ID", "Timepoint"])]
duplicate_mice_df.head()
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Drug Regimen</th>
      <th>Sex</th>
      <th>Age_months</th>
      <th>Weight (g)</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>909</th>
      <td>g989</td>
      <td>Propriva</td>
      <td>Female</td>
      <td>21</td>
      <td>26</td>
      <td>0</td>
      <td>45.000000</td>
      <td>0</td>
    </tr>
    <tr>
      <th>911</th>
      <td>g989</td>
      <td>Propriva</td>
      <td>Female</td>
      <td>21</td>
      <td>26</td>
      <td>5</td>
      <td>47.570392</td>
      <td>0</td>
    </tr>
    <tr>
      <th>913</th>
      <td>g989</td>
      <td>Propriva</td>
      <td>Female</td>
      <td>21</td>
      <td>26</td>
      <td>10</td>
      <td>49.880528</td>
      <td>0</td>
    </tr>
    <tr>
      <th>915</th>
      <td>g989</td>
      <td>Propriva</td>
      <td>Female</td>
      <td>21</td>
      <td>26</td>
      <td>15</td>
      <td>53.442020</td>
      <td>0</td>
    </tr>
    <tr>
      <th>917</th>
      <td>g989</td>
      <td>Propriva</td>
      <td>Female</td>
      <td>21</td>
      <td>26</td>
      <td>20</td>
      <td>54.657650</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>

```python
# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_df = combine_data_df.drop_duplicates(subset=["Mouse ID", "Timepoint"])
html = clean_df.head().to_html()
print(html)
clean_df.head()
```

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Drug Regimen</th>
      <th>Sex</th>
      <th>Age_months</th>
      <th>Weight (g)</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>0</td>
      <td>45.000000</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>5</td>
      <td>38.825898</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>10</td>
      <td>35.014271</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>15</td>
      <td>34.223992</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>k403</td>
      <td>Ramicane</td>
      <td>Male</td>
      <td>21</td>
      <td>16</td>
      <td>20</td>
      <td>32.997729</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>

```python
# Checking the number of mice in the clean DataFrame.
unique_mice_id_count = len(clean_df["Mouse ID"].unique())
print(f"Unique Mice Ids: {unique_mice_id_count}")
```
Unique Mice Ids: 249
  
## Summary Statistics
```python
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
clean_min_df = clean_df[["Drug Regimen", "Tumor Volume (mm3)"]]

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
clean_min_group = clean_min_df.groupby(["Drug Regimen"])

mean_vol = clean_min_group["Tumor Volume (mm3)"].mean().map("{:.2f}".format)
median_vol = clean_min_group["Tumor Volume (mm3)"].median().map("{:.2f}".format)
mode_vol = clean_min_group["Tumor Volume (mm3)"].apply(lambda x: x.mode().iloc[0]).map("{:.2f}".format)
var_vol = clean_min_group["Tumor Volume (mm3)"].var().map("{:.2f}".format)
std_vol = clean_min_group["Tumor Volume (mm3)"].std().map("{:.2f}".format)
sem_vol = clean_min_group["Tumor Volume (mm3)"].sem().map("{:.2f}".format)

# Assemble the resulting series into a single summary dataframe.
summary_df = pd.DataFrame({
                            "Mean": mean_vol,
                            "Median": median_vol,
                            "Mode": mode_vol,
                            "Variance": var_vol,
                            "Standard Deviation": std_vol,
                            "SEM": sem_vol
                         })
summary_df
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mean</th>
      <th>Median</th>
      <th>Mode</th>
      <th>Variance</th>
      <th>Standard Deviation</th>
      <th>SEM</th>
    </tr>
    <tr>
      <th>Drug Regimen</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Capomulin</th>
      <td>40.68</td>
      <td>41.56</td>
      <td>45.00</td>
      <td>24.95</td>
      <td>4.99</td>
      <td>0.33</td>
    </tr>
    <tr>
      <th>Ceftamin</th>
      <td>52.59</td>
      <td>51.78</td>
      <td>45.00</td>
      <td>39.29</td>
      <td>6.27</td>
      <td>0.47</td>
    </tr>
    <tr>
      <th>Infubinol</th>
      <td>52.88</td>
      <td>51.82</td>
      <td>45.00</td>
      <td>43.13</td>
      <td>6.57</td>
      <td>0.49</td>
    </tr>
    <tr>
      <th>Ketapril</th>
      <td>55.24</td>
      <td>53.70</td>
      <td>45.00</td>
      <td>68.55</td>
      <td>8.28</td>
      <td>0.60</td>
    </tr>
    <tr>
      <th>Naftisol</th>
      <td>54.33</td>
      <td>52.51</td>
      <td>45.00</td>
      <td>66.17</td>
      <td>8.13</td>
      <td>0.60</td>
    </tr>
  </tbody>
</table>
<div>

```python
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line
clean_agg_df = clean_min_df.groupby(["Drug Regimen"]).agg({"Tumor Volume (mm3)": ['mean', 'median', lambda x: x.mode().iloc[0], 'var', 'std', 'sem']})
clean_agg_df.rename(columns={"<lambda_0>": "mode"}, inplace=True)
# clean_agg_df[('Tumor Volume (mm3)', 'mean')] = clean_agg_group[('Tumor Volume (mm3)', 'mean')].map("{:.2f}".format)
clean_agg_df
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th colspan="6" halign="left">Tumor Volume (mm3)</th>
    </tr>
    <tr>
      <th></th>
      <th>mean</th>
      <th>median</th>
      <th>mode</th>
      <th>var</th>
      <th>std</th>
      <th>sem</th>
    </tr>
    <tr>
      <th>Drug Regimen</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Capomulin</th>
      <td>40.675741</td>
      <td>41.557809</td>
      <td>45.0</td>
      <td>24.947764</td>
      <td>4.994774</td>
      <td>0.329346</td>
    </tr>
    <tr>
      <th>Ceftamin</th>
      <td>52.591172</td>
      <td>51.776157</td>
      <td>45.0</td>
      <td>39.290177</td>
      <td>6.268188</td>
      <td>0.469821</td>
    </tr>
    <tr>
      <th>Infubinol</th>
      <td>52.884795</td>
      <td>51.820584</td>
      <td>45.0</td>
      <td>43.128684</td>
      <td>6.567243</td>
      <td>0.492236</td>
    </tr>
    <tr>
      <th>Ketapril</th>
      <td>55.235638</td>
      <td>53.698743</td>
      <td>45.0</td>
      <td>68.553577</td>
      <td>8.279709</td>
      <td>0.603860</td>
    </tr>
    <tr>
      <th>Naftisol</th>
      <td>54.331565</td>
      <td>52.509285</td>
      <td>45.0</td>
      <td>66.173479</td>
      <td>8.134708</td>
      <td>0.596466</td>
    </tr>
  </tbody>
</table>
<div>
  

## Bar and Pie Charts
```python
# Generate a bar plot showing the total number of unique mice tested on each drug regimen using pandas.
clean_min_df = clean_df[["Drug Regimen", "Mouse ID"]]
clean_min_group = clean_min_df.groupby(["Drug Regimen"]).count()
clean_min_group.plot(kind="bar")
plt.ylim((0, 235))
plt.title("Total Unique Mice Per Drug Regimen")
plt.xlabel("Drug Regimen")
plt.ylabel("Mice Count")
plt.savefig("./Images/TotalUniqueMiceBar_pandas.png", bbox_inches='tight')
plt.show()
```
![TotalUniqueMiceBar_pandas](Images/TotalUniqueMiceBar_pandas.png)
```python
# Generate a bar plot showing the total number of unique mice tested on each drug regimen using pyplot.
drug_regimen_list = clean_min_group.index.values.tolist()
# print(drug_regimen_list)
plt.bar(drug_regimen_list, clean_min_group["Mouse ID"], label="Mouse ID")
plt.ylim((0, 235))
plt.xticks(rotation="vertical")
plt.title("Total Unique Mice Per Drug Regimen")
plt.xlabel("Drug Regimen")
plt.ylabel("Mice Count")
plt.legend(loc="best")
plt.savefig("./Images/TotalUniqueMiceBar_pyplot.png", bbox_inches='tight')
plt.show()
```
![TotalUniqueMiceBar_pyplot](Images/TotalUniqueMiceBar_pyplot.png)
```python
# Generate a pie plot showing the distribution of female versus male mice using pandas
gender_count = clean_df["Sex"].value_counts()
explode = [0.1, 0]
gender_count.plot(kind="pie", autopct="%1.1f%%", explode=explode, shadow=True)
plt.title("Male versus Female Mouse Distribution")
plt.savefig("./Images/MalevFemaleDist_pandas.png")
plt.show()
```
![MalevFemaleDist_pandas](Images/MalevFemaleDist_pandas.png)
```python
# Generate a pie plot showing the distribution of female versus male mice using pyplot
colors = ["tab:blue", "tab:orange"]
explode = [0.1, 0]
genders = gender_count.index.values.tolist()
plt.pie(gender_count, explode=explode, labels=genders, colors=colors, autopct="%1.1f%%", shadow=True)
plt.title("Male versus Female Mouse Distribution")
plt.savefig("./Images/MalevFemaleDist_pyplot.png")
plt.show()
```
![MalevFemaleDist_pyplot](Images/MalevFemaleDist_pyplot.png)
  
## Quartiles, Outliers and Boxplots
```python
# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
capo_df = clean_df.loc[clean_df["Drug Regimen"] == "Capomulin",:]
rami_df = clean_df.loc[clean_df["Drug Regimen"] == "Ramicane", :]
infu_df = clean_df.loc[clean_df["Drug Regimen"] == "Infubinol", :]
ceft_df = clean_df.loc[clean_df["Drug Regimen"] == "Ceftamin", :]

# Start by getting the last (greatest) timepoint for each mouse
capo_last_group = capo_df.groupby("Mouse ID").max()["Timepoint"]
rami_last_group = rami_df.groupby("Mouse ID").max()["Timepoint"]
infu_last_group = infu_df.groupby("Mouse ID").max()["Timepoint"]
ceft_last_group = ceft_df.groupby("Mouse ID").max()["Timepoint"]
# print(capo_last_group)

# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
capo_last_df = pd.merge(capo_last_group, clean_df, on=("Mouse ID", "Timepoint"), how='left')
rami_last_df = pd.merge(rami_last_group, clean_df, on=("Mouse ID", "Timepoint"), how='left')
infu_last_df = pd.merge(infu_last_group, clean_df, on=("Mouse ID", "Timepoint"), how='left')
ceft_last_df = pd.merge(ceft_last_group, clean_df, on=("Mouse ID", "Timepoint"), how='left')

capo_last_df.head()
```
<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Timepoint</th>
      <th>Drug Regimen</th>
      <th>Sex</th>
      <th>Age_months</th>
      <th>Weight (g)</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>b128</td>
      <td>45</td>
      <td>Capomulin</td>
      <td>Female</td>
      <td>9</td>
      <td>22</td>
      <td>38.982878</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>b742</td>
      <td>45</td>
      <td>Capomulin</td>
      <td>Male</td>
      <td>7</td>
      <td>21</td>
      <td>38.939633</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>f966</td>
      <td>20</td>
      <td>Capomulin</td>
      <td>Male</td>
      <td>16</td>
      <td>17</td>
      <td>30.485985</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>g288</td>
      <td>45</td>
      <td>Capomulin</td>
      <td>Male</td>
      <td>3</td>
      <td>19</td>
      <td>37.074024</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>g316</td>
      <td>45</td>
      <td>Capomulin</td>
      <td>Female</td>
      <td>22</td>
      <td>22</td>
      <td>40.159220</td>
      <td>2</td>
    </tr>
  </tbody>
</table>
</div>

```python
# Put treatments into a list for for loop (and later for plot labels)
labels = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
capo_tumor = capo_last_df["Tumor Volume (mm3)"]
rami_tumor = rami_last_df["Tumor Volume (mm3)"]
infu_tumor = infu_last_df["Tumor Volume (mm3)"]
ceft_tumor = ceft_last_df["Tumor Volume (mm3)"]

# Create list to fill with tumor vol data (for plotting)
tumor_data = [capo_tumor, rami_tumor, infu_tumor, ceft_tumor]

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
    # Locate the rows which contain mice on each drug and get the tumor volumes --> already done
    # add subset --> already done
    # Determine outliers using upper and lower bounds
# Capomulin
capo_quartiles = capo_tumor.quantile([0.25,0.5,0.75])
capo_lowerq = capo_quartiles[0.25]
capo_upperq = capo_quartiles[0.75]
capo_iqr = capo_upperq - capo_lowerq
capo_lower_bound = capo_lowerq - (1.5 * capo_iqr)
capo_upper_bound = capo_upperq + (1.5 * capo_iqr)
print(f"Capomulin outliers could be below {capo_lower_bound:.2f} mm3 or above {capo_upper_bound:.2f} mm3")
# Ramicane
rami_quartiles = rami_tumor.quantile([0.25,0.5,0.75])
rami_lowerq = rami_quartiles[0.25]
rami_upperq = rami_quartiles[0.75]
rami_iqr = rami_upperq - rami_lowerq
rami_lower_bound = rami_lowerq - (1.5 * rami_iqr)
rami_upper_bound = rami_upperq + (1.5 * rami_iqr)
print(f"Ramicane  outliers could be below {rami_lower_bound:.2f} mm3 or above {rami_upper_bound:.2f} mm3")
# Infubinol
infu_quartiles = infu_tumor.quantile([0.25,0.5,0.75])
infu_lowerq = infu_quartiles[0.25]
infu_upperq = infu_quartiles[0.75]
infu_iqr = infu_upperq - infu_lowerq
infu_lower_bound = infu_lowerq - (1.5 * infu_iqr)
infu_upper_bound = infu_upperq + (1.5 * infu_iqr)
print(f"Infubinol outliers could be below {infu_lower_bound:.2f} mm3 or above {infu_upper_bound:.2f} mm3")
# Ceftamin
ceft_quartiles = ceft_tumor.quantile([0.25,0.5,0.75])
ceft_lowerq = ceft_quartiles[0.25]
ceft_upperq = ceft_quartiles[0.75]
ceft_iqr = ceft_upperq - ceft_lowerq
ceft_lower_bound = ceft_lowerq - (1.5 * ceft_iqr)
ceft_upper_bound = ceft_upperq + (1.5 * ceft_iqr)
print(f"Ceftamin  outliers could be below {ceft_lower_bound:.2f} mm3 or above {ceft_upper_bound:.2f} mm3")
```
Capomulin outliers could be below 20.70 mm3 or above 51.83 mm3  
Ramicane  outliers could be below 17.91 mm3 or above 54.31 mm3  
Infubinol outliers could be below 36.83 mm3 or above 82.74 mm3  
Ceftamin  outliers could be below 25.36 mm3 or above 87.67 mm3  
```python
# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
green_diamond = dict(markerfacecolor="g", marker="D")
plt.subplots()
plt.title("Final Tumor Volume Range (mm3) Per Regimen")
plt.xlabel("Regimen")
plt.ylabel("Final Tumor Volumen (mm3)")
plt.boxplot(tumor_data, labels = labels, flierprops=green_diamond)
plt.savefig("./Images/TumorVolumeRange_pyplot.png")
plt.show()
```
![TumorVolumeRange_pyplot](Images/TumorVolumeRange_pyplot.png)
  
## Line and Scatter Plots
```python
# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
mouse_id = "b128"
b128_df = capo_df.loc[capo_df["Mouse ID"] == mouse_id, :]
x_axis = b128_df["Timepoint"]
y_axis = b128_df["Tumor Volume (mm3)"]
plt.plot(x_axis, y_axis)
plt.xlim((0, 45))
plt.ylim((37,46))
plt.title(f"Capomulin Treatment for mouse #{mouse_id}")
plt.xlabel("Time (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.savefig("./Images/TumorVolumeLine_pyplot.png")
plt.show()
```
![TumorVolumeLine_pyplot](Images/TumorVolumeLine_pyplot.png)
```python
# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
capo_average_tumor_vol = capo_df.groupby(["Mouse ID"]).mean()
x_axis = capo_average_tumor_vol["Weight (g)"]
y_axis = capo_average_tumor_vol["Tumor Volume (mm3)"]
plt.scatter(x_axis, y_axis)
plt.xlim((15, 25))
plt.ylim((34, 46))
plt.title("Average Tumor Volume vs. Mouse Weight for the Capomulin Regimen")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.savefig("./Images/TumorVolumeScatter_pyplot.png")
plt.show()
```
![TumorVolumeScatter_pyplot](Images/TumorVolumeScatter_pyplot.png)
  
## Correlation and Regression
```python
# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_axis, y_axis)
print(f"Linear Regression: slope={slope:.2f}, y-intercept={intercept:.2f}, r-value={rvalue:.2f}, p-value={pvalue:.2f}, stderr={stderr:.2f}")
regress_values = x_axis * slope + intercept

plt.scatter(x_axis, y_axis)
plt.xlim((15, 25))
plt.ylim((34, 46))
plt.plot(x_axis, regress_values, "r-")
plt.title("Average Tumor Volume vs. Mouse Weight for the Capomulin Regimen")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.savefig("./Images/TumorVolumeScatterRegress_pyplot.png")
plt.show()
```
Linear Regression: slope=0.95, y-intercept=21.55, r-value=0.84, p-value=0.00, stderr=0.13  
![TumorVolumeScatterRegress_pyplot](Images/TumorVolumeScatterRegress_pyplot.png)