function [Female_data_age,Male_data_age]=count_demographic_information(data)
Age=data.Age;
Gender=data.Gender;
Female_count=0;
Male_count=0;
for i=1:size(Age,1)
    switch char(Gender{i})
        case 'Female'
            Female_count=Female_count+1;
            Female_data_age(Female_count)=Age{i};
        case 'Male'
            Male_count=Male_count+1;
            Male_data_age(Male_count)=Age{i};
    end
end