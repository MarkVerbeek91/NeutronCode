workspace "NeutronCodeWorkspace"
   configurations { "Debug", "Release" }

project "NeutronCode"
   kind "ConsoleApp"
   language "c"
   files { "**.h", "**.c" }
   includedirs{ "include", "include/EnergySpectra", "include/NeutronProduction", "include/ParticleFlux" }
   
   buildoptions { "-std=c99" } 
   
   filter { "configurations:Debug" }
      defines { "DEBUG" }
      flags { "Symbols" }

   filter { "configurations:Release" }
      defines { "NDEBUG" }
      optimize "On"