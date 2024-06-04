const std = @import("std");

pub const File = struct {
    ptr: []align(std.mem.page_size) u8,
    len: u64,
    allocator: std.mem.Allocator,

    pub fn init(fd: std.os.fd_t, allocator: std.mem.Allocator) !File {
        const stats = try std.os.fstat(fd);
        if (stats.size == 0) {
            return error.ZeroFile;
        }
        const ptr = try std.os.mmap(null, @intCast(stats.size), std.os.PROT.READ | std.os.PROT.WRITE, std.os.MAP.SHARED, fd, 0);
        return File{ .ptr = ptr, .len = @intCast(stats.size), .allocator = allocator };
    }

    pub fn deinit(self: *File) void {
        std.os.munmap(self.ptr);
        self.ptr = undefined;
        self.len = 0;
    }
};

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const args = try std.process.argsAlloc(allocator);
    if (args.len != 2) {
        std.debug.print("usage: {s} file\n", .{args[0]});
        std.process.exit(1);
    }
    const path = args[1];

    // mmap input
    var file = try std.fs.cwd().openFile(path, .{ .mode = .read_write });
    defer file.close();
    const len = try file.getEndPos();
    if (len == 0) {
      return;
    }
    const ptr = try std.posix.mmap(null, len, std.posix.PROT.READ | std.posix.PROT.WRITE, .{ .TYPE = .SHARED }, file.handle, 0);
    defer std.posix.munmap(ptr);

    // search for /nix/store
    var i: usize = 0;
    const needle = "/nix/store/";
    while (std.mem.indexOfPos(u8, ptr, i, needle)) |idx| {
        i = idx + needle.len;
        const j = i + 32; // pos of - in a true path
        if (j < len and ptr[j] == '-') {
            std.mem.copyForwards(u8, ptr[i..], "00000000000000000000000000000000");
        }
    }
}
