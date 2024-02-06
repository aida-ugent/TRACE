//destination: "http://backend:8000/backend/:path*",
// use localhost when running without docker

/** @type {import('next').NextConfig} */
const nextConfig = {
  rewrites: async () => {
    return [
      {
        source: "/backend/:path*",
        destination:
          process.env.DOCKER_ENV === "true"
            ? "http://backend:8000/backend/:path*"
            : "http://localhost:8000/backend/:path*"
      },
    ];
  },
};

module.exports = nextConfig;