/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    "./public/index.html",
    "./src/**/*.{js,jsx,ts,tsx}",
  ],
  theme: {
      fontFamily: {
        sans: ["Roboto", "sans-serif"],
      },
      extend: {
        padding: {
          '1/2': '50%',
          full: '100%',
        },
      },
  },
  plugins: [],
}