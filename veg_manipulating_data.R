# Manipulate data for analyses

'%ni%' <- Negate('%in%')
library(dplyr)

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_extracted_data.RData")
veg_manipulated_data <- data.frame(veg_extracted_data)

#trying another way - array
veg_manipulated_data <- veg_manipulated_data %>% mutate(Array_str =
                                                                 case_when(Array < 0 ~ "Axiom", 
                                                                           Array > 0 ~ "BiLEVE")
)

#axiom but numbers
veg_manipulated_data <- veg_manipulated_data %>% mutate(Array_num =
                                                                                case_when(Array < 0 ~ -1, 
                                                                                          Array > 0 ~ 1)
)

#smoking
veg_manipulated_data <- veg_manipulated_data %>% mutate(smoking_standardized =
                                                                                           case_when(Smoking_status == -3 ~ NA, 
                                                                                                     Smoking_status == 0 ~ "Never", 
                                                                                                     Smoking_status == 1 ~ "Previous", 
                                                                                                     Smoking_status == 2 ~ "Current")
)

#physical activity
veg_manipulated_data <- veg_manipulated_data %>% mutate(physical_activity_standardized =
                                                          case_when(Physical_activity == 0 ~ "Low", 
                                                                    Physical_activity == 1 ~ "Moderate", 
                                                                    Physical_activity == 2 ~ "High")
)

#alcohol
veg_manipulated_data <- veg_manipulated_data %>% mutate(alcohol_standardized =
                                                                                                      case_when(Alcohol_status == -3 ~ NA, 
                                                                                                                Alcohol_status == 0 ~ "Never", 
                                                                                                                Alcohol_status == 1 ~ "Previous", 
                                                                                                                Alcohol_status == 2 ~ "Current")
)

#statins
veg_manipulated_data <- veg_manipulated_data %>% mutate(statin_sum = rowSums(across(c(Statin0, Statin1, Statin2, Statin3, Statin4, Statin5, Statin6, Statin7, Statin8, Statin9, Statin10, 
                                                                                                                         Statin11, Statin12, Statin13, Statin14, Statin15, Statin16, Statin17, Statin18, Statin19, Statin20, 
                                                                                                                         Statin21, Statin22, Statin23, Statin24, Statin25, Statin26, Statin27, Statin28, Statin29, Statin30, 
                                                                                                                         Statin31, Statin32, Statin33, Statin34, Statin35, Statin36, Statin37, Statin38, Statin39, Statin40, 
                                                                                                                         Statin41, Statin42, Statin43, Statin44, Statin45, Statin46, Statin47)), na.rm = TRUE))
veg_manipulated_data <- veg_manipulated_data %>% mutate(statin_sum_standardized =
                                                                                                case_when(statin_sum == 0 ~ 0, 
                                                                                                          statin_sum > 0 ~ 1)
)

#assessment centers
veg_manipulated_data <- veg_manipulated_data %>% mutate(ac10003 =
                                                                                            case_when(Assessment_centres == 10003 ~ 1, 
                                                                                                      Assessment_centres != 10003 ~ 0)
)
###
veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11001 =
                                                                            case_when(Assessment_centres == 11001 ~ 1, 
                                                                                      Assessment_centres != 11001 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11002 =
                                                                            case_when(Assessment_centres == 11002 ~ 1, 
                                                                                      Assessment_centres != 11002 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11003 =
                                                                            case_when(Assessment_centres == 11003 ~ 1, 
                                                                                      Assessment_centres != 11003 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11004 =
                                                                            case_when(Assessment_centres == 11004 ~ 1, 
                                                                                      Assessment_centres != 11004 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11005 =
                                                                            case_when(Assessment_centres == 11005 ~ 1, 
                                                                                      Assessment_centres != 11005 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11006 =
                                                                            case_when(Assessment_centres == 11006 ~ 1, 
                                                                                      Assessment_centres != 11006 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11007 =
                                                                            case_when(Assessment_centres == 11007 ~ 1, 
                                                                                      Assessment_centres != 11007 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11008 =
                                                                            case_when(Assessment_centres == 11008 ~ 1, 
                                                                                      Assessment_centres != 11008 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11009 =
                                                                            case_when(Assessment_centres == 11009 ~ 1, 
                                                                                      Assessment_centres != 11009 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11010 =
                                                                            case_when(Assessment_centres == 11010 ~ 1, 
                                                                                      Assessment_centres != 11010 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11011 =
                                                                            case_when(Assessment_centres == 11011 ~ 1, 
                                                                                      Assessment_centres != 11011 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11012 =
                                                                            case_when(Assessment_centres == 11012 ~ 1, 
                                                                                      Assessment_centres != 11012 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11013 =
                                                                            case_when(Assessment_centres == 11013 ~ 1, 
                                                                                      Assessment_centres != 11013 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11014 =
                                                                            case_when(Assessment_centres == 11014 ~ 1, 
                                                                                      Assessment_centres != 11014 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11015 =
                                                                            case_when(Assessment_centres == 11015 ~ 1, 
                                                                                      Assessment_centres != 11015 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11016 =
                                                                            case_when(Assessment_centres == 11016 ~ 1, 
                                                                                      Assessment_centres != 11016 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11017 =
                                                                            case_when(Assessment_centres == 11017 ~ 1, 
                                                                                      Assessment_centres != 11017 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11018 =
                                                                            case_when(Assessment_centres == 11018 ~ 1, 
                                                                                      Assessment_centres != 11018 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11020 =
                                                                            case_when(Assessment_centres == 11020 ~ 1, 
                                                                                      Assessment_centres != 11020 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11021 =
                                                                            case_when(Assessment_centres == 11021 ~ 1, 
                                                                                      Assessment_centres != 11021 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11022 =
                                                                            case_when(Assessment_centres == 11022 ~ 1, 
                                                                                      Assessment_centres != 11022 ~ 0)
)

veg_manipulated_data <- veg_manipulated_data %>% mutate(ac11023 =
                                                                            case_when(Assessment_centres == 11023 ~ 1, 
                                                                                      Assessment_centres != 11023 ~ 0)
)

#save files
veg_manipulated_data <- veg_manipulated_data
save(veg_manipulated_data, file = "veg_manipulated_data.RData")

