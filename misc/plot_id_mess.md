Two trees from Kilwa in two separate plots, each with two stems measured in two consecutive censuses (some columns omitted):

plotcode | year | tag_id | stem_id | multiple |
TKW-001  | 2010 | A00    | 1       | 1        |
TKW-001  | 2012 | A00    | 1       | 1        |
TKW-001  | 2010 | A01    | 2       | 1        |
TKW-001  | 2012 | A01    | 2       | 1        |
TKW-002  | 2010 | A00    | 1       | 1        |
TKW-002  | 2012 | A00    | 1       | 1        |
TKW-002  | 2010 | A01    | 2       | 1        |
TKW-002  | 2012 | A01    | 2       | 1        |

Therefore, to filter to each stem:

```r
s_kilwa <- s %>%
  filter(plotcode %in% kilwa_chosen) %>%
  group_by(plotcode, stem_id) %>%
  filter(year == max(year))
```

And to filter further to each tree:

```r
t_kilwa <- s_kilwa %>%
  group_by(plotcode, multiple) %>%
  summarise(...)
```

In contrast, the setup of trees as recorded in Nhambita:

plotcode | year | tag_id | stem_id | multiple |
MGR-010  | 2006 | 5      | 1841    | 5        |
MGR-010  | 2007 | 5      | 1842    | 5        |
MGR-010  | 2006 | 6      | 1843    | 6        |
MGR-010  | 2007 | 6      | 1844    | 6        |
MGR-011  | 2006 | 5      | 1845    | 5        |
MGR-011  | 2007 | 5      | 1846    | 5        |
MGR-011  | 2006 | 6      | 1847    | 6        |
MGR-011  | 2007 | 6      | 1848    | 6        |

```r
s_nham <- s %>%
  filter(grepl("MGR", plotcode)) %>%
  group_by(plotcode, tag_id) %>%
  filter(year == max(year))
```

And to filter to each tree:

```r
t_nham <- s_nham %>%
  group_by(plotcode, multiple) %>%
  summarise(...)
```

And the same in DRC:

plotcode | year | tag_id | stem_id | multiple |
DKS-002  | 2012	| 10     | 111_1   | NA       |
DKS-002  | 2012	| 10     | 111_2   | NA       |
DKS-002  | 2012 | 11     | 112_1   | NA       |
DKS-002  | 2012 | 11     | 112_2   | NA       |
DKS-003  | 2012	| 10     | 111_1   | NA       |
DKS-003  | 2012	| 10     | 111_2   | NA       |
DKS-003  | 2012 | 11     | 112_1   | NA       |
DKS-003  | 2012 | 11     | 112_2   | NA       |

```r
s_drc <- s %>%
  filter(plotcode %in% c("DKS002", "DKS003")) %>%
  bind_rows(., s_drc_subs) %>%
  group_by(stem_id, plotcode) %>%
  filter(year == max(year)) %>%
  mutate(multiple = gsub("^.*_", "", .$stem_id)) %>%
  mutate(stem_id = gsub("^(.*)[.].*", "\\1", .$stem_id))
```

And to filter to each tree:

```r
t_drc <- s_drc %>%
  group_by(plotcode, stem_id) %>%  
  filter(year == max(year))
```

