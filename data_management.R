source("setup.R")

### GRASP DATA ---------------------------------------------------------------------------------------

UK_recommendations = data.frame(atb=c("Ciprofloxacin","Cefixime OR Ceftriaxone","Ceftriaxone AND Azithromycin"),
                                phase=c("A","B","C"),
                                start=c(2000,2005,2010),
                                end=c(2005,2010,2018))

# Azithromycin --------------------

grasp_azithromycin_mic = 
  read.csv("data/GRASP/azithromycin_mic.csv") %>%
  tbl_df() %>%
  dplyr::select(-prop,-sample_size) %>%
  spread(gender_or,n) %>%
  mutate(HMW=`Heterosexual Male`+Women) %>%
  dplyr::select(year,mic,HMW,MSM) %>%
  gather("gender_or","n",3:4) %>%
  group_by(year,gender_or) %>%
  filter(mic!=512) %>%
  mutate(sample_size=sum(n),
         atb="Azithromycin",
         res=ifelse(mic>1,1,0), # resistance defined as MIC>1.0 mg/L (EUCAST: http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.pdf)
         EUCAST_breakpoint=">1 mg/L",
         micres=ifelse(res==1,EUCAST_breakpoint,paste0(mic," mg/L")),
         micres=factor(micres,levels=c("0.03 mg/L" , "0.06 mg/L" , "0.125 mg/L", "0.25 mg/L" , "0.5 mg/L" ,"1 mg/L",">1 mg/L" ))) %>%
  ungroup()


# Ciprofloxacin --------------------

grasp_ciprofloxacin_mic = 
  read.csv("data/GRASP/ciprofloxacin_mic.csv") %>%
  tbl_df() %>%
  dplyr::select(-prop,-sample_size) %>%
  spread(gender_or,n) %>%
  mutate(HMW=`Heterosexual Male`+Women) %>%
  dplyr::select(year,mic,HMW,MSM) %>%
  gather("gender_or","n",3:4) %>%
  group_by(year,gender_or) %>%
  mutate(sample_size=sum(n),
         atb="Ciprofloxacin",
         res=ifelse(mic>=0.125,1,0), # resistance defined as MIC≥0.125 mg/L (EUCAST: http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.pdf)
         EUCAST_breakpoint=">0.06 mg/L",
         micres=ifelse(res==1,EUCAST_breakpoint,paste0(mic," mg/L")),
         micres=factor(micres,levels=c("0.002 mg/L" ,"0.004 mg/L" ,"0.008 mg/L", "0.015 mg/L", "0.03 mg/L" , "0.06 mg/L" ,">0.06 mg/L"  ))) %>%
  ungroup()

# Cefixime --------------------

grasp_cefixime_mic = 
  read.csv("data/GRASP/cefixime_mic.csv") %>%
  tbl_df() %>%
  dplyr::select(-prop,-sample_size) %>%
  spread(gender_or,n) %>%
  mutate(HMW=`Heterosexual Male`+Women) %>%
  dplyr::select(year,mic,HMW,MSM) %>%
  gather("gender_or","n",3:4) %>%
  group_by(year,gender_or) %>%
  mutate(sample_size=sum(n),
         atb="Cefixime",
         res=ifelse(mic>=0.25,1,0), # resistance defined as MIC≥0.25 mg/L (EUCAST: http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.pdf)
         EUCAST_breakpoint=">0.125 mg/L",
         micres=ifelse(res==1,EUCAST_breakpoint,paste0(mic," mg/L")),
         micres=factor(micres,levels=c("0.002 mg/L",  "0.004 mg/L" , "0.008 mg/L",  "0.015 mg/L" , "0.03 mg/L"  , "0.06 mg/L"  , "0.125 mg/L", ">0.125 mg/L"  ))) %>%
  ungroup()

# Ceftriaxone --------------------

grasp_ceftriaxone_mic = 
  read.csv("data/GRASP/ceftriaxone_mic.csv") %>%
  tbl_df() %>%
  dplyr::select(-prop,-sample_size) %>%
  spread(gender_or,n) %>%
  mutate(HMW=`Heterosexual Male`+Women) %>%
  dplyr::select(year,mic,HMW,MSM) %>%
  gather("gender_or","n",3:4) %>%
  group_by(year,gender_or) %>%
  mutate(sample_size=sum(n),
         atb="Ceftriaxone",
         res=ifelse(mic>=0.25,1,0), # resistance defined as MIC≥0.25 mg/L (EUCAST: http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.pdf)
         EUCAST_breakpoint=">0.125 mg/L",
         micres=ifelse(res==1,EUCAST_breakpoint,paste0(mic," mg/L")),
         micres=factor(micres,levels=c("0.002 mg/L",  "0.004 mg/L" , "0.008 mg/L",  "0.015 mg/L" , "0.03 mg/L"  , "0.06 mg/L"  , "0.125 mg/L", ">0.125 mg/L"  ))) %>%
  ungroup()


# Prescriptions --------------------

grasp_prescriptions =
  read.csv("data/GRASP/grasp_prescriptions.csv") %>%
  tbl_df() %>%
  dplyr::select(-prop,-sample_size) %>%
  spread(gender_or,n) %>%
  mutate(HMW=`Heterosexual Male`+Women) %>%
  dplyr::select(year,treatment,HMW,MSM) %>%
  gather("gender_or","n",3:4) %>%
  group_by(year,gender_or) %>%
  mutate(sample_size=sum(n))

grasp_prescriptions_by_atb =
  grasp_prescriptions %>%
  filter(treatment!="Other",treatment!="Unknown") %>%
  group_by(year,gender_or) %>%
  mutate(sample_size=sum(n),p=n/sample_size,
         azithromycin=ifelse(treatment %in% c("Azithromycin monotherapy","Cefixime and azithromycin","Ceftriaxone and azithromycin"),p,0),
         ciprofloxacin=ifelse(treatment %in% c("Ciprofloxacin monotherapy"),p,0),
         cefixime=ifelse(treatment %in% c("Cefixime and azithromycin","Cefixime Monotherapy"),p,0),   
         ceftriaxone=ifelse(treatment %in% c("Ceftriaxone monotherapy","Ceftriaxone and azithromycin"),p,0)) %>%
  group_by(year,gender_or) %>%
  summarise( Azithromycin=sum(azithromycin),
             Ciprofloxacin=sum(ciprofloxacin),
             Cefixime=sum(cefixime),   
             Ceftriaxone=sum(ceftriaxone)) %>%
  gather("atb","p",3:6) %>%
  ungroup()
         
         
