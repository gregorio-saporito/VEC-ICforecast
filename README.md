# VEC-ICforecast

Inspired by the presentation "Real time forecasting of Covid-19 intensive care units demand" held the 30th of October 2020 as part of the "COVid-19 Empirical Research Webinar", I decided to test the model proposed by Berta P. et al. (2020) to forecast intensive care units demand in Denmark. For more details on the event see <https://ceeds.unimi.it/webinar-cover/>.

The model proposed exploits the cointegrating relationship between the time series of hospitalised patients and intensive care unit (IC) occupancy with a vector error correction (VEC) model. In contrast with an ordinary VAR, the VEC incorporates an error correction term which accounts for any deviations from the long-run equilibrium between the two time series. To read the whole analysis see the R Markdown output which I called :point_right: [Report.pdf](https://github.com/gregorio-saporito/VEC-ICforecast/blob/main/Report.pdf)
