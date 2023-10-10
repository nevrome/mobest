# Parameter estimation for optimal ancestry interpolation

## Using a subset of the variogram to estimate the nugget parameter

samples_projected <- readr::read_csv("docs/data/samples_projected.csv")

### Determining pairwise distances

distances_all <- mobest::calculate_pairwise_distances(
  independent = mobest::create_spatpos(
    id = samples_projected$Sample_ID,
    x = samples_projected$x,
    y = samples_projected$y,
    z = samples_projected$Date_BC_AD_Median
  ),
  dependent = mobest::create_obs(
    C1 = samples_projected$MDS_C1,
    C2 = samples_projected$MDS_C2
  )
)
