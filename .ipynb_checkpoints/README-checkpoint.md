# Connectivity modeling

- need the dev version of `terra`. Download within R via:

```r
remotes::install_github("rspatial/terra",  repos='https://rspatial.r-universe.dev')
remotes::install_github("KennethTM/flowdem")
```

## To dos

- make a git repo
- handle all the hardcoded 10m assumptions (or don't- maybe equal area is fine?)
- right now you are generating figures for every gage (rating curve, inundation gif, floodplain profile). I don't really need to do that for every gage in the US....
- remove 'Hb from final dfs for network model (or rescale it correctly). Right now it is using the Bieger scaling relation, but I've removed that from the analysis and now just use gage-identified Hb. SO maybe upscale that or just completely remove.
- Decide whether to remove the MC analysis or to keep it... I think it makes sense but the comput is gonna go crazzzzy. Maybe only run in some waterhsds For benchmarking for now it's just set to 1 MC using the mean
- in the CT there are ~500 reaches that the model won't run on (thalweg pour points don't overlap the watershed). I think similar to that weird little catchment in the Passumpsic. SHould we care? Not sure this can be easily fixed and they are all probably weird connector reaches
- Scaling below 10m Wb