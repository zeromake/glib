add_rules("mode.debug", "mode.release")

target("glib_gcvt")
    set_kind("$(kind)")
    add_files("src/gcvt.c")
    add_headerfiles("src/gcvt.h", {prefixdir = "glib"})
