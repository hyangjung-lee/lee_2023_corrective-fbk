function [stored_var] = parload_synthetic( filename )

load( filename, 'des' , 'acrR','funInputFurther'); 

stored_var.des = des;
stored_var.acrR = acrR;
stored_var.funInputFurther = funInputFurther;

end
