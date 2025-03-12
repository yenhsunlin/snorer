<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>
<style>
.mono {
    font-family: monospace;
}
</style>


# snorer.params


### *`class`* <span class="mono">snorer.params</span>

This class contains default parameters that will be used in many functions as keyword arguments (`**kwargs`).
Four classes of parameters are predefined and incorporated as attributes in this class. We illustrate them in the followings.



####  <span class="mono">snorer.params.min_distance</span>

Instance.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `d_cut` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Terminating point for $d$. Below this value `snorer.sn_nu_spectrum` will return 0. Default is $3.24\times 10^{-15}$ kpc, approximating 100 km, the size of neutrino sphere. 

> `r_cut` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Terminating $n_\chi$ when $r^\prime <$ `r_cut`, kpc. Below this value `snorer.differential_flux` will return 0. If one needs to incorporate dark matter spike in the central region, `r_cut` cannot be too large. Otherwise, the spike effect will be chopped off before it has any noticeble consequence. Default is $10^{-8}$ kpc.


####  <span class="mono">snorer.params.halo</span>
Instance.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `rhos` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Characteristic density, MeV cm<sup>−3</sup>. Default is 184.

> `rs` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Characteristic length, kpc. Default is 24.4.

> `n` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Index for the halo. Default is 2.
 
#### <span class="mono">snorer.params.spike</span>

Instance.

The following will be used as keyword arguments when `is_spike = True`. When having `is_spike = False`, typing any of the following arguments into function input will result in ValueError.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `mBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Supermassive black hole (SMBH) mass, $M_\odot$. Default is $4.29\times 10^6$.

> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;SMBH age, years. Default is $10^9$.

> `rh` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;SMBH influence radius, kpc. Defaut is 0.002.

> `alpha` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;Profile slope for halo spike, `'3/2'` or `'7/3'`. Default is `'3/2'`.

> `sigv` : *None or float* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM annihilation cross section, in the unit of $10^{-26}$ cm<sup>3</sup> s<sup>−1</sup>. **None** indicates no annihilation and 5.9 as $5.9\times 10^{-26}$.

####  <span class="mono">snorer.params.vegas</span>

Instance.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `nitn` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Number of chains in **vegas** for each integration. Default is 10.

> `neval` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Number of evaluationg points in each chain in **vegas**. Default is 10000.


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

What type or types of parameters are contained in `**kwargs`
of that function will be documented in its docstrings.
Suppose the docstrings says the keywaord arguments containing the values in **min_distance**, **halo** and **spike**, it means all the attributes in the above theree categories can be used as keyword arguments in that function.
We illustrate this in the following example.


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**


Let's take [`snorer.event`](../main/event.md){:target="_blank"} as an example. By looking at its doc page, the necessary `**kwargs` at least contains: `rhos`, `rs`, `n`, `d_cut`, `r_cut`, `nitn` and `neval`. However if one sets `is_spike = True`, then the followings: `mBH`, `tBH`, `rh`, `alpha` and `sigv` are also mandatory. If user didn't specify any of them above, they **snorer** will use default values listed above.

For instance, without any keyword argument specification

```python
>>> import snorer as sn # import snorer
>>> mx,Rs,beta = 1e-2,8.5,0.3
>>> event = sn.event(mx,Rs,beta) # none of the kwargs is specified.
>>> print(event)
1.4301996553521633e-06
```

If we want to specify `rhos` and `mBH`, we can do
```python
>>> event_spike = sn.event(mx,Rs,beta,is_spike=True,rhos=500,mBH=1e7) 
>>> print(event_spike)
3.883407626493457e-06
```

Note that `event_spike` is slightly larger than `event` due to we turn on spike feature. However, the enhancement is not drastic due to supernova is not exactly located at the center as `beta` is not 0. 

If we remove `is_spike`, then the program uses default `is_spike = False`.
This will result in ValueError due to `mBH` only valid when spike feature turns on.