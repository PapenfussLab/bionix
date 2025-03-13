const std = @import("std");

pub const File = struct {
    ptr: []align(std.heap.page_size_min) u8,
    len: u64,
    allocator: std.mem.Allocator,

    pub fn init(fd: std.posix.fd_t, allocator: std.mem.Allocator) !File {
        const stats = try std.posix.fstat(fd);
        if (stats.size == 0) {
            return error.ZeroFile;
        }
        const ptr = try std.posix.mmap(null, @intCast(stats.size), std.posix.PROT.READ | std.posix.PROT.WRITE, .{ .TYPE = .SHARED}, fd, 0);
        return File{ .ptr = ptr, .len = @intCast(stats.size), .allocator = allocator };
    }

    pub fn deinit(self: *File) void {
        std.posix.munmap(self.ptr);
        self.ptr = undefined;
        self.len = 0;
    }
};

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const args = try std.process.argsAlloc(allocator);
    if (args.len == 1) {
        std.process.exit(0);
    } else if (args.len > 2) {
        std.debug.print("usage: {s} file\n", .{args[0]});
        std.process.exit(1);
    }
    const path = args[1];

    // mmap input
    const fd = try std.posix.open(path, .{.ACCMODE = .RDWR, .CREAT = false, .TRUNC = false}, 0);
    var input = File.init(fd, allocator) catch |err| if (err == error.ZeroFile) {
        return;
    } else {
        return err;
    };
    defer input.deinit();

    // search for /nix/store
    var i: usize = 0;
    const needle = "/nix/store/";
    while (std.mem.indexOfPos(u8, input.ptr, i, needle)) |idx| {
        i = idx + needle.len;
        const j = i + 32; // pos of - in a true path
        if (j < input.len and input.ptr[j] == '-') {
            @memcpy(input.ptr[i..], "00000000000000000000000000000000");
        }
    }
}
