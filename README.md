First-year distribution and movement rates of PIT-tagged North East
Atlantic mackerel (Scomber Scombrus) varies with timing, location, and
fish size at release during the spawning season
================

*Contributors: Kotaro Ono<sup>1†</sup>, Sondre Hølleland<sup>1†</sup>,
Aril Slotte<sup>1</sup>,….*

*<sup>1</sup>Institute of Marine Research.*

*<sup>†</sup> Responsible for the code.*

*Correspondance to: <kotaro.uno@hi.no>*

*Nature paper can be found [here](https://www.hi.no).*

This repository contains the necessary code for reproducing results from
the paper *First-year distribution and movement rates of PIT-tagged
North East Atlantic mackerel (Scomber Scombrus) varies with timing,
location, and fish size at release during the spawning season*.

## Data

All use of the data should cite Slotte et. al.(2021).

The data used can be downloaded via API using functions from the
[fishvice/taggart](https://github.com/fishvice/taggart) package. The
package can be installed by

``` r
devtools::install_github("fishvice/taggart", dependencies = FALSE)
```

Due to some no longer maintained dependent packages, we have also copied
the necessary functions to download the data in the
*R/download\_data\_functions.R* file.

The current data in *Data/current.nc* is average monthly surface current
information (Northward and Eastward sea water velocity), from
[Copernicus Marine Service](http://marine.copernicus.eu).

## References

1.  Ono, Kotaro, Sondre Hølleland, Aril Slotte, …, 2021, *First-year
    distribution and movement rates of PIT-tagged North East Atlantic
    mackerel (Scomber Scombrus) varies with timing, location, and fish
    size at release during the spawning season*, Nature, 1-13.
2.  Slotte, Aril … , 2021,…..

## Authors’ github accounts

**Kotaro Ono** - [kotkot](https://github.com/kotkot)

**Sondre Hølleland** - [holleland](https://github.com/holleland)

## License

This project is licensed under the GNU GPLv3 License - see
[LICENSE.md](LICENSE.md) for details. The data is licenced under ….
**Link to data site here**.
