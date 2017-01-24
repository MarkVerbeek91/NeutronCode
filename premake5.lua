workspace "NeutronCodeWorkspace"
   configurations { "Debug", "Release" }

project "NeutronCode"
   kind "ConsoleApp"
   language "c"
   files { "**.h", "**.c" }
   includedirs{ ".", "include", "include/EnergySpectra", "include/FusionReactionRate", "include/NeutronProduction" }
   
   buildoptions { "-std=c99"} 
   links { "m" }
   
   filter { "configurations:Debug" }
      defines { "DEBUG" }
      symbols "On"
	  buildoptions { "-DDEBUG_PARAMETER" }

   filter { "configurations:Release" }
      defines { "NDEBUG" }
      optimize "On"