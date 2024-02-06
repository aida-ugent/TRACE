


export function getHDNeighbors(selectedPoints, kNeighbors, metric) {
    return new Promise((resolve, reject) => {
        if (selectedPoints.length > 0) {
            fetch("/backend/computeHDNeighbors", {
                method: "POST",
                body: JSON.stringify({
                    k: kNeighbors,
                    points: selectedPoints,
                    hd_metric: metric,
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(response => {
                    if (!response.ok) {
                        console.log(`/backend/computeHDNeighbors got HTTP error: Status ${response.status}, Text ${response.statusText}`);
                        resolve(response);
                    } else {
                        response.json().then(data => {
                            let newOpacities = data["binary"]
                            resolve(newOpacities);
                        })
                    }
                })
            }
        else {
            resolve([]);
        }
    })
}